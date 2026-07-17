# fig6_s2_badk.R — Pareto-k diagnostic plot for Study 2 (pre) PSIS-LOO
# Input: models/s2_pareto_k.rds (from the loo diagnostics run; long: model, obs, k)
# Mirrors fig2_badk.R styling. Output: figure/fig6_s2_pareto_k.pdf (appendix)

library(tidyverse)
library(patchwork)

pk <- readRDS("models/s2_pareto_k.rds")

make_k_plot <- function(df, title_label) {
  df <- df |> mutate(kcat = case_when(k > 1.0 ~ "k > 1.0",
                                      k > 0.7 ~ "k > 0.7",
                                      TRUE ~ "OK"))
  n07 <- sum(df$k > 0.7); n10 <- sum(df$k > 1.0); ntot <- nrow(df)
  ann <- sprintf("k > 0.7: %d (%.1f%%)\nk > 1.0: %d (%.1f%%)",
                 n07, 100*n07/ntot, n10, 100*n10/ntot)
  ggplot(df, aes(x = obs, y = k, color = kcat)) +
    geom_point(shape = 3, size = 0.9, stroke = 0.5) +
    geom_hline(yintercept = 0,   linetype = "dotted", color = "grey60",  linewidth = 0.4) +
    geom_hline(yintercept = 0.7, linetype = "dashed", color = "#D9A0A0", linewidth = 0.5) +
    geom_hline(yintercept = 1.0, linetype = "solid",  color = "#B04040", linewidth = 0.5) +
    scale_color_manual(
      values = c("OK" = "#6497B1", "k > 0.7" = "#E6A817", "k > 1.0" = "#CC2936"),
      guide  = "none"
    ) +
    annotate("text", x = Inf, y = Inf, label = ann,
             hjust = 1.05, vjust = 1.1, size = 2.5, fontface = "bold", family = "Menlo") +
    labs(x = "Data point", y = "Pareto shape k", title = title_label) +
    theme_classic(base_size = 10, base_family = "Arial") +
    theme(plot.title = element_text(size = 9, face = "bold"))
}

pa <- make_k_plot(filter(pk, model == "hyp"),   "Hypothesized")
pb <- make_k_plot(filter(pk, model == "alt1"),  "Alternate 1")
pc <- make_k_plot(filter(pk, model == "mine2"), "Hyp + log(t)")
pd <- make_k_plot(filter(pk, model == "mine1"), "Alt1 + log(t)")

fig6 <- (pa | pb) / (pc | pd) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') +
  theme(plot.tag = element_text(size = 10, face = "bold", family = "Arial"))

ggsave("figure/fig6_s2_pareto_k.pdf", fig6,
       width = 6.5, height = 6.0, units = "in", device = cairo_pdf)
cat("wrote figure/fig6_s2_pareto_k.pdf\n")
