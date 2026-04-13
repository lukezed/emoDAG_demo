# run.R

library(cmdstanr)
library(posterior)
library(loo)
library(tidyverse)

# ---- Data prep ---------------------------------------------------------------

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

mod_hyp  <- cmdstan_model("stan/hyp_optimized.stan")
mod_alt1 <- cmdstan_model("stan/alt1_optimized.stan")
mod_m1   <- cmdstan_model("stan/mine1_logtrial.stan")


# ---- Init functions ----------------------------------------------------------

init_hyp <- function() {
  list(
    pop_gamma_g_sd = 0.5, pop_alpha_g_sd = 0.5,
    pop_beta_gp_g_sd = 0.5, pop_beta_e_g_sd = 0.5,
    pop_gamma_p_sd = 0.5, pop_alpha_p_sd = 0.5,
    pop_beta_gp_p_sd = 0.5, pop_beta_e_p_sd = 0.5,
    pop_gamma_e_sd = 0.5, pop_alpha_e_sd = 0.5, pop_beta_gp_e_sd = 0.5,
    sigma_g = 1, sigma_p = 1, sigma_e = 1
  )
}

# Alt1 (NCP): init SDs and sigma positive
init_alt1 <- function() {
  list(
    pop_gamma_g_sd = 0.5, pop_alpha_g_sd = 0.5,
    pop_beta_s_sd  = 0.5, pop_beta_f_sd  = 0.5, pop_beta_e_g_sd = 0.5,
    pop_gamma_p_sd = 0.5, pop_alpha_p_sd = 0.5,
    pop_beta_gp_p_sd = 0.5, pop_beta_e_p_sd = 0.5,
    pop_gamma_e_sd = 0.5, pop_alpha_e_sd = 0.5, pop_beta_gp_e_sd = 0.5,
    sigma_g = 1, sigma_p = 1, sigma_e = 1
  )
}


# ---- Fit & save --------------------------------------------------------------

fit_hyp <- mod_hyp$sample(
  data = sdata, chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 1000,
  adapt_delta = 0.95, max_treedepth = 15, seed = 2026,
  init = init_hyp
)
fit_hyp$save_object("models/fit_hyp.rds")


fit_alt1 <- mod_alt1$sample(
  data = sdata, chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 1000,
  adapt_delta = 0.95, max_treedepth = 20, seed = 2026,
  init = init_alt1
)
fit_alt1$save_object("models/fit_alt1.rds")


fit_m1 <- mod_m1$sample(
  data = sdata, chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 1000,
  adapt_delta = 0.95, max_treedepth = 15, seed = 2026,
  init = init_alt1
)
fit_m1$save_object("models/fit_mine1.rds")

# Though original study set at adapt_delta = 0.99, max_treedepth = 20, 
# the model can work well under 0.95/15.


# ---- Diagnostics -------------------------------------------------------------

fit_hyp$cmdstan_diagnose()
# No problems detected: no divergences, E-BFMI OK, R-hat OK, ESS OK

fit_alt1$cmdstan_diagnose()
# No problems detected: no divergences, E-BFMI OK, R-hat OK, ESS OK

fit_m1$cmdstan_diagnose()
# No problems detected: no divergences, E-BFMI OK, R-hat OK, ESS OK

# ---- WAIC comparison ----------------------------------------------------------

# Original paper used WAIC for model comparison:
#   Hypothesized: 52,934.7   (original; our: 53,484.4 - within SE)
#   Alternate 1:  52,743.3  (original; ours: 53,326.1 — within SE)
#
# All three models show >16% p_waic > 0.4 warnings → WAIC unreliable.

waic_hyp  <- waic(fit_hyp$draws("log_lik", format = "matrix"))
waic_alt1 <- waic(fit_alt1$draws("log_lik", format = "matrix"))
waic_m1   <- waic(fit_m1$draws("log_lik", format = "matrix"))

print(waic_hyp)
print(waic_alt1)
print(waic_m1)


# WAIC values:
#   hyp:   53484.4 (se: 1149.8) (16.7% p_waic > 0.4)
#   alt1:  53326.1 (se: 1176.6)  (17.2% p_waic > 0.4)
#   mine1: 53095.0 (se: 1156.6) (17.4% p_waic > 0.4)


loo_compare(list(hyp = waic_hyp, alt1 = waic_alt1, mine1 = waic_m1))

#       elpd_diff se_diff
# mine1    0.0       0.0
# alt1  -115.5      15.4
# hyp   -194.7      33.0

loo_compare(list(hyp = waic_hyp, alt1 = waic_alt1))

#       elpd_diff se_diff
# alt1     0.0       0.0
# hyp    -79.7      34.0

# ---- LOO comparison ----------------------------------------------------------

loo_hyp  <- fit_hyp$loo(variables = "log_lik")
loo_alt1 <- fit_alt1$loo(variables = "log_lik")
loo_m1   <- fit_m1$loo(variables = "log_lik")

loo_hyp 
loo_alt1
loo_m1

# Pareto k diagnostics:
#   hyp:   k>0.7: 89 (1.3%),  k>1.0: 25 (0.4%)
#   alt1:  k>0.7: 111 (1.6%), k>1.0: 19 (0.3%)
#   mine1: k>0.7: 113 (1.7%), k>1.0: 25 (0.4%)
#
# LOO Pareto k values appear lower than expected because iter_sampling=1000
# (4000 total draws). With higher draws the bad k counts may change.
# Regardless, the autoregressive structure means LOO allows future
# information leakage — LFO-CV is the appropriate method (see run_lfo.R).


loo_compare(list(alt1 = loo_alt1, mine1 = loo_m1, hyp = loo_hyp))
#       elpd_diff se_diff
# mine1    0.0       0.0
# alt1   -99.8      12.4
# hyp   -177.2      33.8


loo_compare(list(alt1 = loo_alt1, hyp = loo_hyp))
#      elpd_diff se_diff
# alt1   0.0       0.0
# hyp  -77.5      33.6

loo_model_weights(list(alt1 = loo_alt1, mine1 = loo_m1, hyp = loo_hyp))

# Method: stacking
# ------
#   weight
# hyp   0.177 
# alt1  0.000 
# mine1 0.823 

loo_model_weights(list(alt1 = loo_alt1, hyp = loo_hyp))

# Method: stacking
# ------
#   weight
# alt1 0.701 
# hyp  0.299 

loo_model_weights(list(alt1 = loo_alt1, mine1 = loo_m1))

# Method: stacking
# ------
#  weight
# alt1  0.000 
# mine1 1.000 

loo_model_weights(list(mine1 = loo_m1, hyp = loo_hyp))

# Method: stacking
# ------
#   weight
# mine1 0.858 
# hyp   0.142 

# ---- Population-level estimates ----------------------------------------------



pop_vars <- function(fit) {
  v <- fit$metadata()$stan_variables
  grep("^(pop_|mu_|sd_|sigma|tau_)", v, value = TRUE) |>
    setdiff(grep("^z", v, value = TRUE))
}

print(fit_hyp$summary(variables = pop_vars(fit_hyp),
                      mean = mean, median = median,
                      q2.5 = ~quantile(.x, 0.025),
                      q97.5 = ~quantile(.x, 0.975)), n = 30)

#                        median   [95% CrI]
# pop_alpha_g           0.858     [0.838, 0.877]   goal autoregression
# pop_beta_gp_g         0.322     [0.299, 0.346]   goal-perf discrepancy -> goal (symmetric)
# pop_beta_e_g          0.078     [0.049, 0.107]   emotion -> goal
# pop_alpha_p           0.596     [0.554, 0.637]   performance autoregression
# pop_beta_gp_p         0.429     [0.382, 0.474]   goal-perf discrepancy -> performance
# pop_beta_e_p         -0.101     [-0.171,-0.037]  emotion -> performance *** SIGNIFICANT ***
# pop_beta_gp_e         0.128     [0.111, 0.145]   discrepancy -> emotion
# sigma_g               0.871     [0.856, 0.887]
# sigma_p               1.530     [1.500, 1.560]
# sigma_e               0.462     [0.454, 0.471]


print(fit_alt1$summary(variables = pop_vars(fit_alt1),
                       mean = mean, median = median,
                       q2.5 = ~quantile(.x, 0.025),
                       q97.5 = ~quantile(.x, 0.975)), n = 30)
#                        median   [95% CrI]            original [95% HDI]
# pop_alpha_g           0.869     [0.850, 0.888]       —
# pop_beta_s            0.273     [0.250, 0.296]       0.273 [0.249, 0.296]  ✓
# pop_beta_f            0.483     [0.427, 0.541]       0.483 [0.427, 0.539]  ✓
# pop_alpha_p           0.596     [0.555, 0.636]       —
# pop_beta_gp_p         0.428     [0.383, 0.475]       0.429 [0.383, 0.476]  ✓
# pop_beta_e_p         -0.104     [-0.174,-0.036]      -0.103 [-0.172,-0.036] ✓
# pop_beta_gp_e         0.128     [0.112, 0.145]       0.128 [0.112, 0.145]  ✓
# pop_beta_e_g          0.067     [0.038, 0.098]       0.067 [0.038, 0.096]  ✓
# sigma_g               0.858     [0.843, 0.875]
# sigma_p               1.530     [1.510, 1.560]
# sigma_e               0.462     [0.454, 0.471]


mine1_vars <- c(pop_vars(fit_m1), "beta_lt_g", "beta_lt_p", "beta_lt_e")
print(fit_m1$summary(variables = mine1_vars,
                     mean = mean, median = median,
                     q2.5 = ~quantile(.x, 0.025),
                     q97.5 = ~quantile(.x, 0.975)), n = 30)
#                        median   [95% CrI]            change from alt1
# pop_alpha_g           0.877     [0.857, 0.896]       ~stable
# pop_beta_s            0.273     [0.249, 0.296]       ~stable
# pop_beta_f            0.491     [0.434, 0.549]       ~stable
# pop_alpha_p           0.500     [0.456, 0.544]       DECREASED (was 0.596)
# pop_beta_gp_p         0.369     [0.322, 0.418]       DECREASED (was 0.428)
# pop_beta_e_p         -0.043     [-0.114, 0.028]      *** NO LONGER SIGNIFICANT *** (was -0.104)
# pop_beta_gp_e         0.128     [0.112, 0.145]       ~stable
# sigma_p               1.510     [1.480, 1.540]       slightly improved
#
# log(trial) fixed effects:
# beta_lt_g            -0.047     [-0.083,-0.011]      goals decrease over time *
# beta_lt_p             0.326     [ 0.262, 0.393]      *** STRONG practice effect ***
# beta_lt_e            -0.017     [-0.037, 0.002]      emotions: marginal (CrI includes 0)

# ---- Key finding -------------------------------------------------------------
# The negative emotion -> performance effect (beta_e_p = -0.104, 95% CrI
# [-0.174, -0.036]) reported in the original paper (both hyp and alt1)
# becomes non-significant (beta_e_p = -0.043, 95% CrI [-0.114, 0.028])
# after controlling for a log(trial) learning curve.
#
# The original authors explained this effect via the "coasting hypothesis"
# (Carver, 2003): positive emotions signal sufficient progress, leading to
# reduced effort. Our results suggest this may instead be a confound —
# emotions trend downward over time while performance trends upward due to
# practice, creating a spurious negative association.
#
# Additionally, beta_lt_e for emotions is marginal (95% CrI barely includes
# zero), while the original 90% CI [-0.033, -0.001] excluded zero. This
# further suggests the emotion time trend is weak.
#
# The goal equation (success/failure asymmetry) and emotion equation are
# robust to this correction.

