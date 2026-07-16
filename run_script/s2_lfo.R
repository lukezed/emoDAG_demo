# s2_lfo.R — exact one-step-ahead LFO-CV for Study 2 (UNBALANCED panel)
#
# Usage: Rscript run_script/s2_lfo.R <group> <model>
#   group: pre | post      model: hyp | alt1 | mine1 | mine2
#   e.g.  Rscript run_script/s2_lfo.R pre mine1
#
# Differences from lfo_fit.R (Study 1, balanced):
#   1. L = 10 (Burkner et al. 2020 ~1/3 initial history), not 5.
#   2. Unbalanced panel: participants have 20-39 sessions. At step t we fit each
#      person's sessions up to min(t, T_i) -- training always includes all ~305
#      people (only their length differs), so the fits stay stable -- and predict
#      session t+1 ONLY for people who have it. The held-out set shrinks with t.
#   3. Stores POINTWISE one-step-ahead log predictive density per (subj, step),
#      NOT just the per-step sum. This is what makes models comparable on an
#      unbalanced panel: two models are always scored on the SAME held-out
#      (subj, session) points, so their paired difference is valid; the
#      observation-level SE (sqrt(n_points) * sd of pointwise diffs, as in
#      loo_compare) handles the varying per-step N correctly. Per-step summed
#      ELPD is NOT comparable across steps (different N) and is reported only
#      as a diagnostic.
#
# Output: models/s2lfo_<group>_<model>.rds  (long: step, subj, elpd_i)
# Runtime: full analysis (MIN_N=1) predicts sessions 11..38 = 28 steps/model,
#   ~30 min/refit -> ~14 h/model, ~2.4 days for all four (local single-refit).

library(cmdstanr)
library(posterior)
library(tidyverse)

# ---- Args --------------------------------------------------------------------
argv  <- commandArgs(trailingOnly = TRUE)
group <- if (length(argv) >= 1) argv[1] else "pre"
mod_name <- if (length(argv) >= 2) argv[2] else "mine1"
stopifnot(group %in% c("pre", "post"),
          mod_name %in% c("hyp", "alt1", "mine1", "mine2"))

stan_file <- c(
  hyp   = "stan/hyp_optimized.stan",
  alt1  = "stan/alt1_optimized.stan",
  mine1 = "stan/mine1_logtrial.stan",
  mine2 = "stan/mine2_logtrial_sym.stan"
)[[mod_name]]
split <- mod_name %in% c("alt1", "mine1")   # split success/failure vs symmetric GPD

L     <- 10L        # Burkner ~1/3 initial history
MIN_N <- 1L         # keep ALL forecasting points (complete analysis, default).
                    # Raising this drops sparse late sessions -- a COMPUTE saver
                    # (a full refit for a handful of people is costly), NOT a
                    # rigour improvement: the observation-level comparison already
                    # handles varying per-step N correctly, so every point is
                    # legitimate information. Set e.g. 100L only if compute-bound.
SEED  <- 2026
cat(sprintf("=== Study 2 LFO: group=%s model=%s (L=%d, min held-out N=%d) ===\n",
            group, mod_name, L, MIN_N))

# ---- Data (Study-1 format; block = re-sequenced session index) ---------------
df <- read_csv(sprintf("data/data2_%s.csv", group), show_col_types = FALSE) |>
  arrange(subj, block)
all_subjs <- sort(unique(df$subj))
N_subj    <- length(all_subjs)
max_block <- max(df$block)

# Cap the last predicted session where held-out N is still meaningful. Late
# sessions have few participants (a full refit each, for ~1-10 people) -- they
# add noise and cost, so drop them. The paired comparison remains valid; we just
# forecast over the well-populated window. Report what is dropped (no silent cap).
heldout_n  <- sapply((L + 1):max_block, function(t) sum(df$block == t))
last_block <- max(((L + 1):max_block)[heldout_n >= MIN_N])  # last session with N >= MIN_N
cat(sprintf("N = %d | max sessions = %d | predict %d..%d (%d steps); dropped %d sparse late sessions (N<%d)\n",
            N_subj, max_block, L + 1, last_block, last_block - L,
            max_block - last_block, MIN_N))

# ---- Helpers -----------------------------------------------------------------
# Training/eval data up to `last_block`. Unbalanced-safe: each person contributes
# only the sessions they have (block <= last_block); everyone has block 1.
make_sdata <- function(last_block) {
  block1 <- df |> filter(block == 1) |> arrange(subj)
  rest   <- df |> filter(block >= 2, block <= last_block) |> arrange(subj, block)
  list(
    N = N_subj, T = nrow(rest),
    id    = as.integer(factor(rest$subj, levels = all_subjs)),
    trial = rest$block,
    g1 = block1$goals, p1 = block1$perfs, e1 = block1$emots,
    g  = rest$goals,   p  = rest$perfs,   e  = rest$emots
  )
}

make_init <- function() {
  base <- list(
    pop_gamma_g_sd = 0.5, pop_alpha_g_sd = 0.5, pop_beta_e_g_sd = 0.5,
    pop_gamma_p_sd = 0.5, pop_alpha_p_sd = 0.5,
    pop_beta_gp_p_sd = 0.5, pop_beta_e_p_sd = 0.5,
    pop_gamma_e_sd = 0.5, pop_alpha_e_sd = 0.5, pop_beta_gp_e_sd = 0.5,
    sigma_g = 1, sigma_p = 1, sigma_e = 1
  )
  if (split) { base$pop_beta_s_sd <- 0.5; base$pop_beta_f_sd <- 0.5 }
  else       { base$pop_beta_gp_g_sd <- 0.5 }
  base
}

log_mean_exp <- function(x) { m <- max(x); m + log(mean(exp(x - m))) }

# ---- Compile -----------------------------------------------------------------
mod <- cmdstan_model(stan_file)

# ---- LFO loop: fit on 1..t-1, predict session t -----------------------------
pw <- tibble(step = integer(), subj = integer(), elpd_i = double())  # pointwise store

for (t in (L + 1):last_block) {
  t0 <- Sys.time()

  fit <- mod$sample(
    data = make_sdata(t - 1),
    chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 1000,
    adapt_delta = 0.95, max_treedepth = 15,
    seed = SEED, init = make_init, refresh = 0
  )
  if (sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent) > 0)
    cat(sprintf("  [block %d] WARNING: divergences\n", t))

  # one-step-ahead: run gq on data through block t, score the block-t rows only
  sdata_t <- make_sdata(t)
  gq <- mod$generate_quantities(fit, data = sdata_t, seed = SEED)
  idx <- which(sdata_t$trial == t)              # people who HAVE session t
  subj_t <- all_subjs[sdata_t$id[idx]]
  ll   <- gq$draws("log_lik", format = "matrix")[, idx, drop = FALSE]  # draws x n_t
  elpd_i <- apply(ll, 2, log_mean_exp)          # per-person one-step-ahead ELPD

  pw <- bind_rows(pw, tibble(step = t, subj = subj_t, elpd_i = elpd_i))
  saveRDS(pw, sprintf("models/s2lfo_%s_%s.rds", group, mod_name))  # crash-safe

  cat(sprintf("  block %2d: n=%3d  sumELPD=%.1f  [%.1f min]\n",
              t, length(idx), sum(elpd_i),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

cat(sprintf("\nDONE. %s/%s: total ELPD_LFO = %.2f over %d points -> models/s2lfo_%s_%s.rds\n",
            group, mod_name, sum(pw$elpd_i), nrow(pw), group, mod_name))
