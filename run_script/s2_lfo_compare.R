# s2_lfo_compare.R — compare Study 2 LFO across models at the OBSERVATION level
#
# Usage: Rscript run_script/s2_lfo_compare.R <group> [model ...]
#   group: pre | post   (default pre); models default to all four
#
# Why observation-level: on the unbalanced panel two models are scored on the
# same held-out (step, subj) points, so we join on those keys and compute the
# paired difference per point. Total ΔELPD = sum of pointwise diffs; SE =
# sqrt(n_points) * sd(pointwise diffs) -- the loo_compare estimator, valid for
# varying per-step N. (Per-step summed ELPD is NOT comparable across steps.)

library(tidyverse)

argv  <- commandArgs(trailingOnly = TRUE)
group <- if (length(argv) >= 1) argv[1] else "pre"
mods  <- if (length(argv) >= 2) argv[-1] else c("hyp", "alt1", "mine1", "mine2")

# long pointwise: one row per (model, step, subj)
pw <- map_dfr(mods, function(m) {
  f <- sprintf("models/s2lfo_%s_%s.rds", group, m)
  if (!file.exists(f)) { message("missing: ", f); return(NULL) }
  readRDS(f) |> mutate(model = m)
})
present <- unique(pw$model)

cat(sprintf("=== Study 2 (%s) LFO totals (higher = better) ===\n", group))
pw |> group_by(model) |>
  summarise(elpd = sum(elpd_i), n_points = n()) |>
  arrange(desc(elpd)) |> as.data.frame() |> print(row.names = FALSE)

# paired, observation-level: join two models on (step, subj)
cmp <- function(a, b) {
  if (!all(c(a, b) %in% present)) return(invisible())
  j <- inner_join(filter(pw, model == a), filter(pw, model == b),
                  by = c("step", "subj"), suffix = c(".a", ".b"))
  d  <- j$elpd_i.a - j$elpd_i.b
  se <- sd(d) * sqrt(length(d))
  cat(sprintf("%-6s - %-6s:  ΔELPD = %7.2f,  SE = %6.2f  (%.1f SE)  n=%d\n",
              a, b, sum(d), se, sum(d) / se, length(d)))
}

cat("\n=== 2x2 pairwise (paired on matched held-out points) ===\n")
cmp("mine1", "mine2")   # asymmetry, both +logT  <- does asymmetry hold under LFO here?
cmp("mine1", "alt1")    # log-trial, split
cmp("mine2", "hyp")     # log-trial, symmetric
cmp("alt1",  "hyp")     # asymmetry, no logT

# ---- RESULTS (fill in after s2_lfo.R has run for all models) -----------------
# group=pre:
#   (pending run)
# group=post:
#   (pending run)
