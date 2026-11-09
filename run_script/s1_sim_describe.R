# s1_sim_describe.R — compare simulated datasets with the real Study 1 data
# (Appendix D). Regenerates the datasets of s1_sim_check.R (same generator,
# same seeds, no model fits) and writes models/sim_describe.csv.
# Usage: Rscript run_script/s1_sim_describe.R
suppressMessages({library(posterior); library(tidyverse)})
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

sim_once <- function(seed, beta_true) {
  set.seed(seed); out <- vector("list", N)
  for (i in 1:N) {
    g <- p <- e <- numeric(20)
    g[1] <- b1$goals[i]; p[1] <- b1$perfs[i]; e[1] <- b1$emots[i]
    for (t in 2:20) {
      d <- p[t-1]-g[t-1]; lt <- log(t)
      g[t] <- clamp(rnorm(1, ind$gamma_g[i]+ind$alpha_g[i]*g[t-1]+ind$beta_gp_s[i]*max(d,0)+
                      ind$beta_gp_f[i]*min(d,0)+ind$beta_e_g[i]*e[t-1]+pop["beta_lt_g"]*lt, pop["sigma_g"]), 0, 60)
      p[t] <- clamp(rnorm(1, ind$gamma_p[i]+ind$alpha_p[i]*p[t-1]+ind$beta_gp_p[i]*(g[t]-p[t-1])+
                      beta_true*e[t-1]+pop["beta_lt_p"]*lt, pop["sigma_p"]), 0, 50)
      e[t] <- clamp(rnorm(1, ind$gamma_e[i]+ind$alpha_e[i]*e[t-1]+ind$beta_gp_e[i]*(p[t]-g[t])+
                      pop["beta_lt_e"]*lt, pop["sigma_e"]), -3, 3)
    }
    out[[i]] <- tibble(subj=b1$subj[i], block=1:20, goals=g, perfs=p, emots=e)
  }
  bind_rows(out)
}

describe <- function(d) {
  w <- d |> group_by(subj) |>
    mutate(e_lag=lag(emots), lt=log(block)) |>
    mutate(across(c(perfs, emots, lt), ~ .x - mean(.x), .names="c_{.col}")) |>
    filter(!is.na(e_lag)) |> mutate(c_elag=e_lag-mean(e_lag)) |> ungroup()
  m <- function(v, b) mean(d[[v]][d$block==b])
  tibble(P1=m("perfs",1), P10=m("perfs",10), P20=m("perfs",20),
         G1=m("goals",1), G10=m("goals",10), G20=m("goals",20),
         E1=m("emots",1), E10=m("emots",10), E20=m("emots",20),
         sdP=sd(d$perfs), sdE=sd(d$emots),
         r_lt_P=cor(w$c_lt, w$c_perfs), r_lt_E=cor(w$c_lt, w$c_emots),
         r_Elag_P=cor(w$c_elag, w$c_perfs),
         edgeE=mean(abs(d$emots)>=3))
}

real <- describe(df) |> mutate(set="real", rep=0L)
sims <- map_dfr(c(0, -0.10), function(b) map_dfr(1:100, function(r)
  describe(sim_once(2026+r, b)) |> mutate(set=paste0("beta_", b), rep=r)))
res <- bind_rows(real, sims)
write_csv(res, "models/sim_describe.csv")
