# compare_lfo.R — Compare LFO results across models
# Run after all three lfo jobs finish

library(tidyverse)

results <- map(c("hyp", "alt1", "mine1"), ~ readRDS(sprintf("models/lfo_%s.rds", .x)))
wide <- bind_rows(results) |> pivot_wider(names_from = model, values_from = elpd)

cat("Per-block ELPD:\n")
print(wide)

cat("\n--- Totals ---\n")
totals <- bind_rows(results) |> group_by(model) |> summarise(elpd_lfo = sum(elpd))
print(totals)

for (ref in c("hyp", "alt1")) {
  d <- wide$mine1 - wide[[ref]]
  cat(sprintf("\nmine1 - %s:  ΔELPD = %.2f,  SE = %.2f\n",
              ref, sum(d), sd(d) * sqrt(length(d))))
}
d2 <- wide$alt1 - wide$hyp
cat(sprintf("\nalt1 - hyp:   ΔELPD = %.2f,  SE = %.2f\n",
            sum(d2), sd(d2) * sqrt(length(d2))))

write_csv(wide, "models/lfo_comparison.csv")


# ---- LFO-CV Results ---------------------------------------------------------
#
# Per-block ELPD:
#   block    hyp    alt1   mine1
#       6 -1473  -1471  -1467
#       7 -1433  -1437  -1436
#       8 -1610  -1605  -1615
#       9 -1385  -1381  -1381
#      10 -2039  -2061  -2051
#      11 -1427  -1429  -1424
#      12 -1369  -1360  -1354
#      13 -1367  -1365  -1364
#      14 -1279  -1277  -1271
#      15 -1331  -1327  -1323
#      16 -1371  -1369  -1366
#      17 -1335  -1335  -1330
#      18 -1373  -1366  -1361
#      19 -1377  -1369  -1363
#      20 -1752  -1733  -1738
#
# Totals:
#   hyp:   -21922
#   alt1:  -21883
#   mine1: -21844
#
# Pairwise comparisons:
#   mine1 - hyp:   ΔELPD = 78.12,  SE = 28.92  (2.7 SE → strong evidence)
#   mine1 - alt1:  ΔELPD = 39.54,  SE = 19.20  (2.1 SE → sufficient evidence)
#   alt1  - hyp:   ΔELPD = 38.57,  SE = 33.80  (1.1 SE → insufficient evidence)
#
# Key observations:
#   1. mine1 is the best model under LFO-CV, confirming LOO/WAIC results.
#   2. The advantage of mine1 over alt1 shrinks from LOO (99.8) to LFO (39.54)
#      because LOO inflates mine1's advantage via future information leakage —
#      the log(trial) trend benefits most from seeing future data.
#   3. alt1 vs hyp is no longer clearly distinguishable (1.1 SE), suggesting
#      that controlling for the learning curve matters more than splitting
#      success/failure in the goal equation.
# ---- End of results ---------------------------------------------------------
