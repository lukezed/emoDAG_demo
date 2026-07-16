# lfo_compare.R — Compare LFO results across the 2x2 models
# Run after all four lfo jobs finish (hyp, alt1, mine1, mine2)

library(tidyverse)

results <- map(c("hyp", "alt1", "mine1", "mine2"), ~ readRDS(sprintf("models/lfo_%s.rds", .x)))
wide <- bind_rows(results) |> pivot_wider(names_from = model, values_from = elpd)

cat("Per-block ELPD:\n"); print(wide)

cat("\n--- Totals (higher = better) ---\n")
totals <- bind_rows(results) |> group_by(model) |>
  summarise(elpd_lfo = sum(elpd)) |> arrange(desc(elpd_lfo))
print(totals)

# pairwise ΔELPD with paired SE = sd(blockwise diffs) * sqrt(n_blocks)
cmp <- function(a, b) {
  d  <- wide[[a]] - wide[[b]]
  se <- sd(d) * sqrt(length(d))
  cat(sprintf("%-6s - %-6s:  ΔELPD = %7.2f,  SE = %6.2f  (%.1f SE)\n",
              a, b, sum(d), se, sum(d) / se))
}
cat("\n--- 2x2 pairwise: asymmetry (split/sym) x log-trial (with/without) ---\n")
cmp("mine1", "mine2")   # split vs symmetric, both +logT  <- is asymmetry supported once practice modelled?
cmp("mine1", "alt1")    # +logT vs no-logT (split)
cmp("mine2", "hyp")     # +logT vs no-logT (symmetric)
cmp("alt1",  "hyp")     # split vs symmetric, both no-logT (the original contrast)
cmp("mine1", "hyp")     # best vs worst

write_csv(wide, "models/lfo_comparison.csv")


# ---- LFO-CV Results (Study 1, 15 one-step-ahead steps, blocks 6-20) ---------
#
# Totals (higher = better):
#   mine1 (split + logT)  -21843.78   <- best
#   alt1  (split, no logT) -21883.32
#   mine2 (sym   + logT)  -21905.88
#   hyp   (sym,  no logT)  -21921.90   <- worst
#
# 2x2 pairwise (ΔELPD, paired SE):
#   mine1 - mine2:  62.10  SE 20.76  (3.0 SE)  asymmetry, both +logT   -> supported
#   mine1 - alt1 :  39.54  SE 19.20  (2.1 SE)  log-trial, split        -> supported
#   mine2 - hyp  :  16.02  SE 21.41  (0.7 SE)  log-trial, symmetric    -> inconclusive
#   alt1  - hyp  :  38.57  SE 33.80  (1.1 SE)  asymmetry, no logT      -> inconclusive
#   mine1 - hyp  :  78.12  SE 28.92  (2.7 SE)
#
# Key observations:
#   1. mine1 (split + logT) is best under exact LFO, matching LOO/WAIC ordering.
#   2. mine1 - alt1 shrinks from LOO (99.8) to LFO (39.5): LOO inflated the
#      log-trial advantage via future-information leakage.
#   3. THE ASYMMETRY REFRAME. Without a practice term, split vs symmetric
#      (alt1 - hyp) is inconclusive under LFO (1.1 SE) -- consistent with the
#      original Alt1-over-Hyp preference being weak. BUT once a practice term is
#      added, split beats symmetric credibly (mine1 - mine2 = 62.1, 3.0 SE):
#      controlling for practice sharpens the goal-equation contrast and the
#      success/failure asymmetry IS supported. Paper stance: CONFIRM the
#      asymmetry (Lin's Alt1 goal equation), OVERTURN the emotion-to-performance
#      conclusion (see s1_fit.R: beta_E->P non-credible in mine1 & mine2).
# ---- End of results ---------------------------------------------------------
