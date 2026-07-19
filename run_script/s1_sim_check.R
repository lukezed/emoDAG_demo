# s1_sim_check.R — generative check for the confound mechanism (Appendix)
# Simulate from the fitted LogT (mine1) posterior with beta_E->P := 0, fit
# Alt1 (no practice term): does a spurious ~ -0.10 emerge?
# Usage: Rscript run_script/s1_sim_check.R [n_reps]
#
# Simulated trajectories are truncated to generous physical bounds
# (g in [0,60], p in [0,50], e in [-3,3]). The fitted linear system is
# dynamically explosive for a few participants' median parameters
# (unbounded simulation reaches |g| ~ 4e4), while the real task is
# bounded by block time and the rating scale; truncation mimics those
# bounds and keeps the fit geometry sane. Disclose in the appendix.
suppressMessages({library(cmdstanr); library(posterior); library(tidyverse)})
args <- commandArgs(trailingOnly=TRUE); n_reps <- if (length(args)) as.integer(args[1]) else 3L

df <- read_csv("data/data.csv", show_col_types=FALSE) |> arrange(subj, block)
b1 <- df |> filter(block==1) |> arrange(subj); N <- nrow(b1)

fit <- readRDS("models/fit_mine1.rds")
med <- function(v) apply(fit$draws(v, format="matrix"), 2, median)
ind <- sapply(c("gamma_g","alpha_g","beta_gp_s","beta_gp_f","beta_e_g",
                "gamma_p","alpha_p","beta_gp_p","gamma_e","alpha_e","beta_gp_e"),
              med, simplify=FALSE)
pop <- sapply(c("beta_lt_g","beta_lt_p","beta_lt_e","sigma_g","sigma_p","sigma_e"),
              function(v) median(fit$draws(v, format="matrix")))
rm(fit); gc()

clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

sim_once <- function(seed) {
  set.seed(seed); out <- vector("list", N)
  for (i in 1:N) {
    g <- p <- e <- numeric(20)
    g[1] <- b1$goals[i]; p[1] <- b1$perfs[i]; e[1] <- b1$emots[i]
    for (t in 2:20) {
      d <- p[t-1]-g[t-1]; lt <- log(t)
      g[t] <- clamp(rnorm(1, ind$gamma_g[i]+ind$alpha_g[i]*g[t-1]+ind$beta_gp_s[i]*max(d,0)+
                      ind$beta_gp_f[i]*min(d,0)+ind$beta_e_g[i]*e[t-1]+pop["beta_lt_g"]*lt, pop["sigma_g"]), 0, 60)
      p[t] <- clamp(rnorm(1, ind$gamma_p[i]+ind$alpha_p[i]*p[t-1]+ind$beta_gp_p[i]*(g[t]-p[t-1])+
                      0*e[t-1]+pop["beta_lt_p"]*lt, pop["sigma_p"]), 0, 50)   # beta_E->P = 0
      e[t] <- clamp(rnorm(1, ind$gamma_e[i]+ind$alpha_e[i]*e[t-1]+ind$beta_gp_e[i]*(p[t]-g[t])+
                      pop["beta_lt_e"]*lt, pop["sigma_e"]), -3, 3)
    }
    out[[i]] <- tibble(subj=b1$subj[i], block=1:20, goals=g, perfs=p, emots=e)
  }
  bind_rows(out)
}

mod <- cmdstan_model("stan/alt1_optimized.stan")
init_a <- function() list(pop_gamma_g_sd=0.5,pop_alpha_g_sd=0.5,pop_beta_s_sd=0.5,
  pop_beta_f_sd=0.5,pop_beta_e_g_sd=0.5,pop_gamma_p_sd=0.5,pop_alpha_p_sd=0.5,
  pop_beta_gp_p_sd=0.5,pop_beta_e_p_sd=0.5,pop_gamma_e_sd=0.5,pop_alpha_e_sd=0.5,
  pop_beta_gp_e_sd=0.5,sigma_g=1,sigma_p=1,sigma_e=1)

res <- map_dfr(seq_len(n_reps), function(r) {
  sd_ <- sim_once(2026+r); rest <- sd_ |> filter(block>=2)
  sdata <- list(N=N, T=nrow(rest), id=as.integer(factor(rest$subj)), trial=rest$block,
    g1=sd_$goals[sd_$block==1], p1=sd_$perfs[sd_$block==1], e1=sd_$emots[sd_$block==1],
    g=rest$goals, p=rest$perfs, e=rest$emots)
  f <- mod$sample(data=sdata, chains=4, parallel_chains=4, iter_warmup=1000,
    iter_sampling=1000, adapt_delta=0.95, max_treedepth=15, seed=2026, init=init_a, refresh=0)
  q <- quantile(f$draws("pop_beta_e_p", format="matrix"), c(.025,.5,.975))
  cat(sprintf("rep %d: beta_E->P = %.3f [%.3f, %.3f]\n", r, q[2], q[1], q[3]))
  tibble(rep=r, median=q[2], lo=q[1], hi=q[3])
})
write_csv(res, "models/sim_check_results.csv")
cat("DONE. True value 0; original spurious estimate -0.104.\n")

# ---- RESULTS (10 reps, seeds 2027-2036; fits seed 2026, 0.95/15) -------------
# medians: -0.082 -0.061 -0.016 -0.028 -0.052 -0.096 -0.028 +0.014 -0.043 -0.070
# negative in 9/10 (median -0.048); credibly negative in 4/10
# (reps 1, 2, 6, 10; largest -0.096 [-0.154, -0.037]).
# -> with a TRUE null emotion->perf path, omitting the practice term biases
#    the estimate negative by ~ -0.05 on average, matching the real-data
#    attenuation when logT is added (-0.104 -> -0.043). Variability across
#    reps reflects the weak emotion drift that opens the back door
#    (beta_lt_e ~ -0.017). NOTE: unbounded simulation from the median
#    parameters is dynamically explosive for a few participants (|g| ~ 4e4);
#    trajectories are truncated to physical bounds, touched for ~0.3% of
#    g/p values and ~4% of e values.
