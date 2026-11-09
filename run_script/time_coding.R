# time_coding.R — robustness of beta_E->P to the coding of practice time
# Usage: Rscript run_script/time_coding.R <s1|s2pre> <log|linear|dummies>
#   log reproduces mine1 (sanity check); dummies = one level per block from
#   block 3 on (block 2 is the reference), Study 2 blocks >= 30 pooled
#   because few participants reach them.
suppressMessages({library(cmdstanr); library(posterior); library(tidyverse)})
args <- commandArgs(trailingOnly = TRUE); study <- args[1]; coding <- args[2]
f <- c(s1 = "data/data.csv", s2pre = "data/data2_pre.csv")[[study]]
df <- read_csv(f, show_col_types = FALSE) |> arrange(subj, block)
block1 <- df |> filter(block == 1) |> arrange(subj)
rest   <- df |> filter(block >= 2) |> arrange(subj, block)
t <- rest$block
X <- switch(coding,
  log     = cbind(log(t)),
  linear  = cbind(t),
  dummies = { b <- factor(pmin(t, 30)); model.matrix(~ b)[, -1, drop = FALSE] })
sdata <- list(N = n_distinct(df$subj), T = nrow(rest),
  id = as.integer(factor(rest$subj)), trial = t,
  g1 = block1$goals, p1 = block1$perfs, e1 = block1$emots,
  g = rest$goals, p = rest$perfs, e = rest$emots, K = ncol(X), X = unname(X))
init <- function() list(
  pop_gamma_g_sd = 0.5, pop_alpha_g_sd = 0.5, pop_beta_s_sd = 0.5,
  pop_beta_f_sd = 0.5, pop_beta_e_g_sd = 0.5, pop_gamma_p_sd = 0.5,
  pop_alpha_p_sd = 0.5, pop_beta_gp_p_sd = 0.5, pop_beta_e_p_sd = 0.5,
  pop_gamma_e_sd = 0.5, pop_alpha_e_sd = 0.5, pop_beta_gp_e_sd = 0.5,
  sigma_g = 1, sigma_p = 1, sigma_e = 1)
mod <- cmdstan_model("stan/mine1_timeX.stan")
fit <- mod$sample(data = sdata, chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.95,
  max_treedepth = 15, seed = 2026, init = init, refresh = 0)
q <- quantile(fit$draws("pop_beta_e_p", format = "matrix"), c(.025, .5, .975))
s <- fit$summary(c("pop_beta_e_p", "b_p", "b_e"))
print(s, n = Inf)
cat(sprintf("%s %s: beta_E->P = %.3f [%.3f, %.3f], max rhat %.3f, divergences %d\n",
  study, coding, q[2], q[1], q[3], max(fit$summary()$rhat, na.rm = TRUE),
  sum(fit$diagnostic_summary()$num_divergent)))
write_csv(tibble(study, coding, median = q[2], lo = q[1], hi = q[3]),
  "models/time_coding_results.csv", append = TRUE)  # header pre-created

# ---- RESULTS (2026-10-29): beta_E->P median [95% CrI] ------------------------
#                 Study 1                  Study 2 (pre)
# log t           -0.040 [-0.114, 0.029]   (mine1: 0.020)
# linear t        -0.047 [-0.118, 0.022]    0.021 [0.004, 0.037]
# block dummies   -0.043 [-0.114, 0.028]    0.020 [0.004, 0.037]
# log t reproduces mine1 (-0.043 [-0.114, 0.028]) up to MC error. All fits:
# max rhat <= 1.019, 0 divergences. Conclusion is insensitive to the coding:
# Study 1 not credible under every coding, Study 2 credibly positive.
