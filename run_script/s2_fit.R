# s2_fit.R — Study 2 (pre-/post-learning emotions), 2x2 model space
#
# Reuses the Study 1 Stan models unchanged (Study 2 data is in Study-1 format).
#   alt1  = split GPD,     no logT   <- baseline, recovers Lin's beta_e_p = +0.018
#   mine1 = split GPD,     + logT
#   hyp   = symmetric GPD, no logT
#   mine2 = symmetric GPD, + logT
#
# Usage: Rscript run_script/s2_fit.R <pre|post> [model ...]
#   default models: all four (alt1 first). Fits saved as fit_s2<group>_<name>.rds
#
# ---------------------------------------------------------------------------
# RESULTS (population beta_E->P and beta_logt_P; 4 chains, 1000+1000, seed 2026)
#
# PRE group (N=303, T=8553) -- primary; pipeline reproduces Lin's pre baseline:
#   sigma 0.343/0.388 exact in the pre-shift validation run; aligned run 0.341/0.385, alpha_p ~0.70, beta_gp_p ~0.63 == Lin
#
#   model          beta_E->P [95% CrI]        beta_logt_P [95% CrI]
#   alt1 (base)    +0.017 [ 0.001, 0.033]     --
#   mine1 (+logT)  +0.020 [ 0.004, 0.037]     0.025 [0.013, 0.036]
#   hyp  (base)    +0.017 [ 0.001, 0.034]     --
#   mine2 (+logT)  +0.020 [ 0.003, 0.036]     0.025 [0.013, 0.037]
#
#   -> the positive effect is ROBUST to log(session): stays credibly positive
#      (and robust to split/symmetric GPD). Standardised beta_E->P ~ +0.028.
#      Contrast Study 1, where the negative effect is NOT robust (attenuates to
#      non-credible). Practice gain on performance is comparable across studies
#      (S2 ~+8-9 attempted Qs, ~+10%; S1 ~+12%). Do not compare raw beta shifts
#      across studies (different scales); report robustness, not shift size.
#
# POST group (N=295) -- negative control (Lin post baseline beta_E->P ~ -0.004 ns):
#   alt1 (base)    -0.007 [-0.019, 0.004]     --
#   mine1 (+logT)  -0.006 [-0.018, 0.005]     0.023 [0.012, 0.035]
#   -> null baseline stays null under log(session); practice term manufactures
#      no spurious effect.
#
# Convergence: all population means R-hat < 1.01; two individual-level SD
# hyperparameters (pop_gamma_p_sd, pop_beta_e_p_sd) reach ~1.01 in some fits,
# not affecting population inference.
# ---------------------------------------------------------------------------

library(cmdstanr)
library(posterior)
library(loo)
library(tidyverse)

# ---- Args: group + optional model filter -------------------------------------
argv  <- commandArgs(trailingOnly = TRUE)
group <- if (length(argv) >= 1) argv[1] else "pre"
stopifnot(group %in% c("pre", "post"))
sel <- argv[-1]

# ---- Data (Study-1 format from s2_prep.R) --------------------------------
df <- read_csv(sprintf("data/data2_%s.csv", group), show_col_types = FALSE) |>
  arrange(subj, block)
block1 <- df |> filter(block == 1) |> arrange(subj)
rest   <- df |> filter(block >= 2) |> arrange(subj, block)

sdata <- list(
  N = n_distinct(df$subj), T = nrow(rest),
  id = as.integer(factor(rest$subj)), trial = rest$block,
  g1 = block1$goals, p1 = block1$perfs, e1 = block1$emots,
  g = rest$goals, p = rest$perfs, e = rest$emots
)
cat(sprintf("Study 2 (%s): N = %d, T = %d\n", group, sdata$N, sdata$T))

# ---- Model table -------------------------------------------------------------
models <- tribble(
  ~name,   ~stan,                            ~split,
  "alt1",  "stan/alt1_optimized.stan",       TRUE,
  "mine1", "stan/mine1_logtrial.stan",       TRUE,
  "hyp",   "stan/hyp_optimized.stan",        FALSE,
  "mine2", "stan/mine2_logtrial_sym.stan",   FALSE
)
if (length(sel)) models <- filter(models, name %in% sel)

# init SDs positive; success/failure vs symmetric goal SD depending on model
init_for <- function(split) function() {
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

pop_vars <- function(fit) {
  v <- fit$metadata()$stan_variables
  grep("^(pop_|sigma|beta_lt)", v, value = TRUE) |> setdiff(grep("^z", v, value = TRUE))
}

# ---- Fit each model ----------------------------------------------------------
for (i in seq_len(nrow(models))) {
  m <- models[i, ]
  cat(sprintf("\n===== %s (%s) =====\n", m$name, m$stan))
  fit <- cmdstan_model(m$stan)$sample(
    data = sdata, chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 1000,
    adapt_delta = 0.95, max_treedepth = 15, seed = 2026,
    init = init_for(m$split)
  )
  fit$save_object(sprintf("models/fit_s2%s_%s.rds", group, m$name))
  fit$cmdstan_diagnose()

  vars <- pop_vars(fit)
  est <- fit$summary(variables = vars, median = median,
                     q2.5 = ~quantile(.x, 0.025), q97.5 = ~quantile(.x, 0.975))
  # headline: emotion -> performance (the thesis parameter)
  ep <- est |> filter(variable == "pop_beta_e_p")
  cat(sprintf("  beta_E->P = %.3f [%.3f, %.3f]   (Lin Study 2 pre: +0.018)\n",
              ep$median, ep$q2.5, ep$q97.5))
  if ("beta_lt_p" %in% vars) {
    lp <- est |> filter(variable == "beta_lt_p")
    cat(sprintf("  beta_logt_P = %.3f [%.3f, %.3f]\n", lp$median, lp$q2.5, lp$q97.5))
  }
  print(est, n = 30)
  print(waic(fit$draws("log_lik", format = "matrix")))
}
cat("\nDONE. Fits saved to models/fit_s2_*.rds\n")
