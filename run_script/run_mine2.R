# run_mine2.R — fit Hyp + log(trial) (symmetric GPD), same settings as run.R
#
# Completes the 2x2 model space: {symmetric, split GPD} x {no logtrial, logtrial}
#   hyp   = symmetric,  no logtrial
#   alt1  = split,      no logtrial
#   mine1 = split,      logtrial
#   mine2 = symmetric,  logtrial   <- this script
#
# Usage: Rscript run_script/run_mine2.R

library(cmdstanr)
library(posterior)
library(loo)
library(tidyverse)

# ---- Data prep (identical to run.R) ------------------------------------------

df <- read_csv("data/data.csv", show_col_types = FALSE) |>
  arrange(subj, block)

block1 <- df |> filter(block == 1) |> arrange(subj)
rest   <- df |> filter(block >= 2) |> arrange(subj, block)

sdata <- list(
  N     = n_distinct(df$subj),
  T     = nrow(rest),
  id    = as.integer(factor(rest$subj)),
  trial = rest$block,
  g1    = block1$goals,
  p1    = block1$perfs,
  e1    = block1$emots,
  g     = rest$goals,
  p     = rest$perfs,
  e     = rest$emots
)

cat(sprintf("N = %d, T = %d\n", sdata$N, sdata$T))

# ---- Compile -----------------------------------------------------------------

mod_m2 <- cmdstan_model("stan/mine2_logtrial_sym.stan")

# ---- Init (hyp structure: 11 individual params, symmetric) -------------------

init_m2 <- function() {
  list(
    pop_gamma_g_sd = 0.5, pop_alpha_g_sd = 0.5,
    pop_beta_gp_g_sd = 0.5, pop_beta_e_g_sd = 0.5,
    pop_gamma_p_sd = 0.5, pop_alpha_p_sd = 0.5,
    pop_beta_gp_p_sd = 0.5, pop_beta_e_p_sd = 0.5,
    pop_gamma_e_sd = 0.5, pop_alpha_e_sd = 0.5, pop_beta_gp_e_sd = 0.5,
    sigma_g = 1, sigma_p = 1, sigma_e = 1
  )
}

# ---- Fit & save (same MCMC settings as run.R) --------------------------------

fit_m2 <- mod_m2$sample(
  data = sdata, chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 1000,
  adapt_delta = 0.95, max_treedepth = 15, seed = 2026,
  init = init_m2
)
fit_m2$save_object("models/fit_mine2.rds")

# ---- Diagnostics -------------------------------------------------------------

fit_m2$cmdstan_diagnose()

# ---- WAIC / LOO (self) -------------------------------------------------------

waic_m2 <- waic(fit_m2$draws("log_lik", format = "matrix"))
loo_m2  <- fit_m2$loo(variables = "log_lik")
print(waic_m2)
print(loo_m2)

# ---- Population estimates + log(trial) effects -------------------------------

pop_vars <- function(fit) {
  v <- fit$metadata()$stan_variables
  grep("^(pop_|mu_|sd_|sigma|tau_)", v, value = TRUE) |>
    setdiff(grep("^z", v, value = TRUE))
}

m2_vars <- c(pop_vars(fit_m2), "beta_lt_g", "beta_lt_p", "beta_lt_e")
print(fit_m2$summary(variables = m2_vars,
                     mean = mean, median = median,
                     q2.5 = ~quantile(.x, 0.025),
                     q97.5 = ~quantile(.x, 0.975)), n = 30)

# ---- Comparison against existing fits ----------------------------------------
# Reuses saved fits from run.R. The key contrast: does the symmetric-GPD model,
# once given a log(trial) term, match or beat the split-GPD alt1/mine1? If so,
# the success/failure asymmetry is not supported with or without practice.

cmp <- tryCatch({
  loo_hyp  <- readRDS("models/fit_hyp.rds")$loo(variables = "log_lik")
  loo_alt1 <- readRDS("models/fit_alt1.rds")$loo(variables = "log_lik")
  loo_m1   <- readRDS("models/fit_mine1.rds")$loo(variables = "log_lik")

  cat("\n=== LOO compare: all four models ===\n")
  print(loo_compare(list(hyp = loo_hyp, alt1 = loo_alt1,
                         mine1 = loo_m1, mine2 = loo_m2)))

  cat("\n=== mine2 (sym+logT) vs mine1 (split+logT) ===\n")
  print(loo_compare(list(mine1 = loo_m1, mine2 = loo_m2)))

  cat("\n=== stacking weights, all four ===\n")
  print(loo_model_weights(list(hyp = loo_hyp, alt1 = loo_alt1,
                               mine1 = loo_m1, mine2 = loo_m2)))
  TRUE
}, error = function(e) { cat("Comparison skipped:", conditionMessage(e), "\n"); FALSE })

cat("\nDONE. Fit saved to models/fit_mine2.rds\n")
cat("beta_e_p and beta_lt_p are the parameters to read first.\n")
