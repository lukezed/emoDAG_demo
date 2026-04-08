# run_lfo.R — Exact Leave-Future-Out CV (1-step-ahead)
#
# Usage: Rscript run_lfo.R <model_name>
#   model_name: hyp, alt1, or mine1

library(cmdstanr)
library(posterior)
library(tidyverse)

# ---- Which model? ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
mod_name <- if (length(args) >= 1) args[1] else "mine1"
stopifnot(mod_name %in% c("hyp", "alt1", "mine1"))

stan_file <- c(
  hyp   = "stan/hyp_optimized.stan",
  alt1  = "stan/alt1_optimized.stan",
  mine1 = "stan/mine1_logtrial.stan"
)[[mod_name]]

cat(sprintf("=== LFO-CV for model: %s ===\n", mod_name))

# ---- Settings ----------------------------------------------------------------

L    <- 5L
SEED <- 2026

# ---- Data prep ---------------------------------------------------------------

df <- read_csv("data/data.csv", show_col_types = FALSE) |>
  arrange(subj, block)

all_subjs <- sort(unique(df$subj))
N_subj    <- length(all_subjs)
max_block <- max(df$block)

cat(sprintf("Subjects: %d | Blocks: %d | L: %d | Steps: %d\n",
            N_subj, max_block, L, max_block - L))

# ---- Helpers -----------------------------------------------------------------

make_sdata <- function(last_block) {
  block1 <- df |> filter(block == 1) |> arrange(subj)
  rest   <- df |> filter(block >= 2, block <= last_block) |> arrange(subj, block)
  list(
    N = N_subj, 
    T = nrow(rest),
    id    = as.integer(factor(rest$subj, levels = all_subjs)),
    trial = rest$block,
    g1 = block1$goals, 
    p1 = block1$perfs, 
    e1 = block1$emots,
    g  = rest$goals,   
    p  = rest$perfs,   
    e  = rest$emots
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
  if (mod_name == "hyp") {
    base$pop_beta_gp_g_sd <- 0.5
  } else {
    base$pop_beta_s_sd <- 0.5
    base$pop_beta_f_sd <- 0.5
  }
  base
}

log_mean_exp <- function(x) {
  mx <- max(x)
  mx + log(mean(exp(x - mx)))
}

# ---- Compile -----------------------------------------------------------------

mod <- cmdstan_model(stan_file)

# ---- LFO loop ----------------------------------------------------------------

elpd_blocks <- numeric(0)
block_ids   <- integer(0)

for (t in (L + 1):max_block) {
  cat(sprintf("\n--- Block %d: fit on 1..%d, predict %d ---\n", t, t - 1, t))
  t0 <- Sys.time()
  
  # 1. Fit on blocks 1..t-1
  fit <- mod$sample(
    data = make_sdata(t - 1),
    chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 1000,
    adapt_delta = 0.95, max_treedepth = 15,
    seed = SEED, init = make_init, refresh = 100
  )
  
  diag <- fit$diagnostic_summary(quiet = TRUE)
  if (sum(diag$num_divergent) > 0)
    cat(sprintf("  WARNING: %d divergences\n", sum(diag$num_divergent)))
  
  # 2. Run generated quantities on data 1..t using posterior from 1..t-1
  sdata_t <- make_sdata(t)
  gq <- mod$generate_quantities(fit, data = sdata_t, seed = SEED)
  
  # 3. Find which columns correspond to block t
  #    Data is sorted by (subj, block), so block t rows are scattered.
  #    sdata_t$trial gives the block number for each row.
  block_t_idx <- which(sdata_t$trial == t)
  stopifnot(length(block_t_idx) == N_subj)
  
  # 4. Extract log_lik for those rows only
  ll_all <- gq$draws("log_lik", format = "matrix")   # S x T_rows
  ll_block_t <- ll_all[, block_t_idx]                # S x N_subj
  
  # 5. ELPD = sum over subjects of log(mean(exp(log_lik)))
  elpd_t <- sum(apply(ll_block_t, 2, log_mean_exp))
  
  elpd_blocks <- c(elpd_blocks, elpd_t)
  block_ids   <- c(block_ids, t)
  
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("  ELPD = %.2f  [%.1f min]\n", elpd_t, elapsed))
  
  # Save after every step (crash protection)
  saveRDS(
    tibble(model = mod_name, block = block_ids, elpd = elpd_blocks),
    sprintf("models/lfo_%s.rds", mod_name)
  )
}

cat(sprintf("\n%s  ELPD_LFO = %.2f\n", mod_name, sum(elpd_blocks)))
cat("Done.\n")

