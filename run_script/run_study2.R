# run_study2.R — Study 2 (pre-learning emotions), 2x2 model space
#
# Reuses the Study 1 Stan models unchanged (Study 2 data is in Study-1 format).
#   alt1  = split GPD,     no logT   <- baseline, must recover Lin's beta_e_p = +0.018
#   mine1 = split GPD,     + logT
#   hyp   = symmetric GPD, no logT
#   mine2 = symmetric GPD, + logT
#
# Usage: Rscript run_script/run_study2.R <pre|post> [model ...]
#   default models: all four (alt1 first). Fits saved as fit_s2<group>_<name>.rds

library(cmdstanr)
library(posterior)
library(loo)
library(tidyverse)

# ---- Args: group + optional model filter -------------------------------------
argv  <- commandArgs(trailingOnly = TRUE)
group <- if (length(argv) >= 1) argv[1] else "pre"
stopifnot(group %in% c("pre", "post"))
sel <- argv[-1]

# ---- Data (Study-1 format from prep_study2.R) --------------------------------
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
