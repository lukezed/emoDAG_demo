# logo.R — PSIS leave-one-group-out (person) CV for the 2x2 model space
#
# Target: a new learner (cold start given block-1 values), complementing LFO
# (next block of observed learners). Person effects are integrated out
# analytically: every equation is linear in its person-level coefficients
# given the observed regressors, coefficients are independent Gaussians, and
# residuals are independent across equations, so per person and equation
#   y ~ N(Z mu + beta_lt * log t,  Z D Z' + sigma^2 I),  D = diag(sd^2).
# Computed by Woodbury from the existing fits (Riha et al., 2026 with the
# exact Gaussian marginal). Same procedure as pointwise_ic_leakage arm 4.
#
# Persons with Pareto k > 0.7 are repaired by exact refits without that
# person (the protocol's branch for at most 10 flagged per candidate).
#
# Usage: Rscript run_script/logo.R <s1|s2pre> [refit]
#   refit: also run the exact refits (~32 min each in S1, ~50 in S2); each
#   result is appended to models/logo_<study>_refit.csv so a rerun resumes.
suppressMessages({library(posterior); library(loo); library(tidyverse); library(parallel)})
argv <- commandArgs(trailingOnly = TRUE); study <- argv[1]; do_refit <- "refit" %in% argv
stopifnot(study %in% c("s1", "s2pre"))

df <- read_csv(if (study == "s1") "data/data.csv" else "data/data2_pre.csv",
               show_col_types = FALSE) |> arrange(subj, block)
b1 <- df |> filter(block == 1) |> arrange(subj)
rest <- df |> filter(block >= 2) |> arrange(subj, block)
persons <- split(rest, factor(rest$subj, levels = b1$subj))

# Woodbury log density of y ~ N(m, sigma^2 I + Z diag(d2) Z'); r = y - m
lmvn_wb <- function(r, Z, ZtZ, d2, s2) {
  n <- length(r); k <- ncol(Z)
  L <- chol(ZtZ + diag(s2 / d2, k))
  b <- backsolve(L, crossprod(Z, r), transpose = TRUE)
  quad <- (sum(r^2) - sum(b^2)) / s2
  logdet <- (n - k) * log(s2) + sum(log(d2)) + 2 * sum(log(diag(L)))
  -0.5 * (n * log(2 * pi) + logdet + quad)
}

# one runnable check: Woodbury == direct dense MVN
local({
  set.seed(1); Z <- cbind(1, rnorm(6)); d2 <- c(.5, .2); s2 <- .3; r <- rnorm(6)
  S <- s2 * diag(6) + Z %*% diag(d2) %*% t(Z); L <- chol(S)
  direct <- -0.5 * (6 * log(2 * pi) + 2 * sum(log(diag(L))) +
                    sum(backsolve(L, r, transpose = TRUE)^2))
  stopifnot(abs(direct - lmvn_wb(r, Z, crossprod(Z), d2, s2)) < 1e-10)
})

# design per person: regressors use observed lags, block 1 as initial condition
design <- function(d, b, split) {
  pg <- c(b$goals, head(d$goals, -1)); pp <- c(b$perfs, head(d$perfs, -1))
  pe <- c(b$emots, head(d$emots, -1)); disc <- pp - pg
  Zg <- if (split) cbind(1, pg, pmax(disc, 0), pmin(disc, 0), pe) else cbind(1, pg, disc, pe)
  list(lt = log(d$block),
       g = list(y = d$goals, Z = Zg),
       p = list(y = d$perfs, Z = cbind(1, pp, d$goals - pp, pe)),
       e = list(y = d$emots, Z = cbind(1, pe, d$perfs - d$goals)))
}

coefs <- function(split) list(
  g = c("gamma_g", "alpha_g", if (split) c("beta_s", "beta_f") else "beta_gp_g", "beta_e_g"),
  p = c("gamma_p", "alpha_p", "beta_gp_p", "beta_e_p"),
  e = c("gamma_e", "alpha_e", "beta_gp_e"))

models <- c(hyp = FALSE, alt1 = TRUE, mine2 = FALSE, mine1 = TRUE)   # value = split GPD
labels <- c(hyp = "Hyp", alt1 = "Alt1", mine2 = "Hyp+log t", mine1 = "Alt1+log t")
has_lt <- c(hyp = FALSE, alt1 = FALSE, mine2 = TRUE, mine1 = TRUE)

draws_of <- function(fit, m) {
  cf <- coefs(models[[m]])
  vars <- c(paste0("pop_", unlist(cf)), paste0("pop_", unlist(cf), "_sd"),
            "sigma_g", "sigma_p", "sigma_e", if (has_lt[[m]]) c("beta_lt_g", "beta_lt_p", "beta_lt_e"))
  da <- fit$draws(vars)
  list(chain = rep(seq_len(nchains(da)), each = niterations(da)),
       dr = as.matrix(as.data.frame(as_draws_matrix(da))))  # plain matrix: draws_matrix rows break diag()
}

# S-vector of log p(person i | draw s), person effects integrated out
ll_person <- function(i, dr, m) {
  cf <- coefs(models[[m]]); S <- nrow(dr); out <- numeric(S)
  ds <- design(persons[[i]], b1[i, ], models[[m]])
  for (eq in c("g", "p", "e")) {
    Z <- ds[[eq]]$Z; ZtZ <- crossprod(Z)
    mu <- dr[, paste0("pop_", cf[[eq]]), drop = FALSE]
    d2 <- dr[, paste0("pop_", cf[[eq]], "_sd"), drop = FALSE]^2
    s2 <- dr[, paste0("sigma_", eq)]^2
    blt <- if (has_lt[[m]]) dr[, paste0("beta_lt_", eq)] else numeric(S)
    R <- ds[[eq]]$y - Z %*% t(mu) - outer(ds$lt, blt)   # n x S residuals
    for (s in seq_len(S)) out[s] <- out[s] + lmvn_wb(R[, s], Z, ZtZ, d2[s, ], s2[s])
  }
  out
}

# exact refit without person j, same sampler settings as s1_fit.R / s2_fit.R
refit_without <- function(j, m) {
  suppressMessages(library(cmdstanr))
  keep <- rest |> filter(subj != b1$subj[j]); kb <- b1[-j, ]
  sdata <- list(N = nrow(kb), T = nrow(keep), id = as.integer(factor(keep$subj, levels = kb$subj)),
                trial = keep$block, g1 = kb$goals, p1 = kb$perfs, e1 = kb$emots,
                g = keep$goals, p = keep$perfs, e = keep$emots)
  init <- function() {
    v <- c(paste0("pop_", unlist(coefs(models[[m]])), "_sd"))
    c(setNames(as.list(rep(0.5, length(v))), v), list(sigma_g = 1, sigma_p = 1, sigma_e = 1))
  }
  # exe_file: reuse the compiled binary; two processes compiling at once race (all chains die)
  fit <- cmdstan_model(exe_file = sub("\\.stan$", "", stan_file[[m]]))$sample(
    data = sdata, chains = 4, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 1000,
    adapt_delta = 0.95, max_treedepth = 15, seed = 2026, init = init, refresh = 0)
  ll <- ll_person(j, draws_of(fit, m)$dr, m)
  log(mean(exp(ll - max(ll)))) + max(ll)
}
stan_file <- c(hyp = "stan/hyp_optimized.stan", alt1 = "stan/alt1_optimized.stan",
               mine2 = "stan/mine2_logtrial_sym.stan", mine1 = "stan/mine1_logtrial.stan")
refit_csv <- sprintf("models/logo_%s_refit.csv", study)

logo <- list()
for (m in names(models)) {
  path <- sprintf("models/fit_%s%s.rds", if (study == "s1") "" else "s2pre_", m)
  cat("loading", path, "\n")
  fit <- readRDS(path); dd <- draws_of(fit, m); rm(fit); gc(verbose = FALSE)
  ll <- do.call(cbind, mclapply(seq_along(persons), ll_person, dr = dd$dr, m = m, mc.cores = 8))
  r_eff <- relative_eff(exp(ll - apply(ll, 2, max)), chain_id = dd$chain)
  logo[[m]] <- loo(ll, r_eff = r_eff)
  k <- logo[[m]]$diagnostics$pareto_k
  cat(sprintf("%-11s elpd_LOGO = %.1f (SE %.1f), k>0.7: %d of %d persons, max k = %.2f\n",
              labels[[m]], logo[[m]]$estimates["elpd_loo", 1], logo[[m]]$estimates["elpd_loo", 2],
              sum(k > 0.7), length(k), max(k)))
  if (do_refit) for (j in which(k > 0.7)) {
    done <- if (file.exists(refit_csv)) read_csv(refit_csv, show_col_types = FALSE) else tibble(model = character(), person = integer())
    if (any(done$model == m & done$person == j)) next
    e <- refit_without(j, m)
    cat(sprintf("  refit %s person %d: PSIS %.2f -> exact %.2f\n", m, j, logo[[m]]$pointwise[j, "elpd_loo"], e))
    write_csv(tibble(model = m, person = j, psis = logo[[m]]$pointwise[j, "elpd_loo"], exact = e),
              refit_csv, append = file.exists(refit_csv))
  }
  if (file.exists(refit_csv)) {
    rf <- read_csv(refit_csv, show_col_types = FALSE) |> filter(model == m)
    logo[[m]]$pointwise[rf$person, "elpd_loo"] <- rf$exact
    cat(sprintf("  %d of %d flagged persons replaced by exact refits\n", nrow(rf), sum(k > 0.7)))
  }
}
saveRDS(logo, sprintf("models/logo_%s.rds", study))

# paired contrasts, SE over persons (the LOGO unit)
pair <- function(a, b) {
  d <- logo[[b]]$pointwise[, "elpd_loo"] - logo[[a]]$pointwise[, "elpd_loo"]
  cat(sprintf("%-11s -> %-11s  diff = %7.1f (SE %5.1f)  z = %5.2f\n",
              labels[[a]], labels[[b]], sum(d), sqrt(length(d)) * sd(d), sum(d) / (sqrt(length(d)) * sd(d))))
}
for (m in names(logo)) cat(sprintf("%-11s elpd_LOGO (after repairs) = %.1f\n", labels[[m]], sum(logo[[m]]$pointwise[, "elpd_loo"])))
cat("\nPractice term:\n"); pair("hyp", "mine2"); pair("alt1", "mine1")
cat("Goal asymmetry:\n");  pair("hyp", "alt1"); pair("mine2", "mine1")

# ---- RESULTS (PSIS-LOGO, flagged persons replaced by exact refits) ----------
# Flagged (k > 0.7): S1 4/4/3/6, S2 2/4/3/3 of 357 / 303 persons (Hyp/Alt1/
# Hyp+log t/Alt1+log t); all replaced. PSIS was close for most; it was
# optimistic for one outlier per study (S1 person 92 by 30-35 elpd, S2
# person 286 by 2-8), which lowers every model about equally.
#
#                           Study 1                 Study 2 (pre)
# Practice term
#   Hyp  -> Hyp+log t       48.7 (13.3) z 3.67      36.3 (10.7) z 3.40
#   Alt1 -> Alt1+log t      42.6 (13.8) z 3.08      34.1 (10.3) z 3.32
# Goal asymmetry
#   Hyp  -> Alt1            88.1 (29.8) z 2.96      32.7 (14.1) z 2.32
#   Hyp+log t -> Alt1+log t 82.0 (32.6) z 2.51      30.5 (16.5) z 1.85
# -> practice term credible in all four; asymmetry credible in three, marginal
#    in S2 once log t is in. LOGO targets a new learner, LFO the next block.
