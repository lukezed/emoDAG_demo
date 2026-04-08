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
  iter_warmup = 1000, iter_sampling = 10000,
  adapt_delta = 0.99, max_treedepth = 20, seed = 2026,
  init = init_hyp
)
fit_hyp$save_object("models/fit_hyp.rds")


fit_alt1 <- mod_alt1$sample(
  data = sdata, chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 10000,
  adapt_delta = 0.99, max_treedepth = 20, seed = 2026,
  init = init_alt1
)
fit_alt1$save_object("models/fit_alt1.rds")



fit_m1 <- mod_m1$sample(
  data = sdata, chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 10000,
  adapt_delta = 0.99, max_treedepth = 20, seed = 2026,
  init = init_alt1
)
fit_m1$save_object("models/fit_mine1.rds")

# Though I set at adapt_delta = 0.99, max_treedepth = 20, the model can work 
# well under 0.95/15.


# ---- Diagnostics -------------------------------------------------------------

fit_hyp$cmdstan_diagnose()
# No problems detected: no divergences, E-BFMI OK, R-hat OK, ESS OK

fit_alt1$cmdstan_diagnose()
# No problems detected: no divergences, E-BFMI OK, R-hat OK, ESS OK

fit_m1$cmdstan_diagnose()
# No problems detected: no divergences, E-BFMI OK, R-hat OK, ESS OK


# ---- LOO comparison ----------------------------------------------------------

loo_hyp  <- fit_hyp$loo(variables = "log_lik")
loo_alt1 <- fit_alt1$loo(variables = "log_lik")
loo_m1   <- fit_m1$loo(variables = "log_lik")
loo_compare(list(alt1 = loo_alt1, mine1 = loo_m1, hyp = loo_hyp))
#       elpd_diff se_diff
# mine1   0.0       0.0
# alt1  -99.8      12.4
#
# mine1 substantially outperforms alt1 (delta_LOOIC ~ 200)

# ---- Population-level estimates ----------------------------------------------

pop_vars <- function(fit) {
  v <- fit$metadata()$stan_variables
  grep("^(pop_|mu_|sd_|sigma|tau_)", v, value = TRUE) |>
    setdiff(grep("^z", v, value = TRUE))
}

print(fit_alt1$summary(variables = pop_vars(fit_alt1)), n = 30)
#                       mean     [q5,    q95]
# pop_alpha_g          0.869     [0.853, 0.885]   goal autoregression
# pop_beta_s           0.273     [0.253, 0.292]   goal adjustment after success
# pop_beta_f           0.483     [0.436, 0.532]   goal adjustment after failure
# pop_alpha_p          0.596     [0.563, 0.631]   performance autoregression
# pop_beta_gp_p        0.428     [0.390, 0.467]   goal-perf discrepancy -> performance
# pop_beta_e_p        -0.104     [-0.161,-0.048]  emotion -> performance *** SIGNIFICANT ***
# pop_beta_gp_e        0.129     [0.115, 0.143]   discrepancy -> emotion
# sigma_g              0.858                       residual SD (goals)
# sigma_p              1.530                       residual SD (performance)
# sigma_e              0.462                       residual SD (emotion)

mine1_vars <- c(pop_vars(fit_m1), "beta_lt_g", "beta_lt_p", "beta_lt_e")
print(fit_m1$summary(variables = mine1_vars), n = 30)
#                       mean     [q5,    q95]       change from alt1
# pop_alpha_g          0.877     [0.861, 0.893]     ~stable
# pop_beta_s           0.273     [0.253, 0.293]     ~stable
# pop_beta_f           0.491     [0.442, 0.538]     ~stable
# pop_alpha_p          0.500     [0.464, 0.537]     DECREASED (was 0.596)
# pop_beta_gp_p        0.369     [0.330, 0.410]     DECREASED (was 0.428)
# pop_beta_e_p        -0.043     [-0.101,+0.017]    *** NO LONGER SIGNIFICANT *** (was -0.104)
# pop_beta_gp_e        0.128     [0.114, 0.142]     ~stable
# sigma_p              1.510                         slightly improved
#
# log(trial) fixed effects:
# beta_lt_g           -0.047     [-0.077,-0.016]    goals decrease slightly over time
# beta_lt_p            0.326     [ 0.272, 0.381]    *** STRONG practice effect ***
# beta_lt_e           -0.017     [-0.033,-0.001]    emotions decrease slightly over time
#
# ---- Key finding -------------------------------------------------------------
# The negative emotion -> performance effect (beta_e_p = -0.104) reported
# in the original paper becomes non-significant (beta_e_p = -0.043, CI
# crosses zero) after controlling for a log(trial) learning curve.
#
# The original authors explained this effect via the "coasting hypothesis"
# (Carver, 2003): positive emotions signal sufficient progress, leading to
# reduced effort. Our results suggest this may instead be a confound —
# emotions trend downward over time while performance trends upward due to
# practice, creating a spurious negative association.
#
# The goal equation (success/failure asymmetry) and emotion equation are
# robust to this correction.


