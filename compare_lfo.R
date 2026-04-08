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
