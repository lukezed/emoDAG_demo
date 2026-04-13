library(tidyverse)
library(extrafont)

lfo <- read_csv("models/lfo_comparison.csv", show_col_types = FALSE)

delta_df <- tibble(
  block = lfo$block,
  `LogT − Alt1`  = lfo$mine1 - lfo$alt1,
  `LogT − Hyp`   = lfo$mine1 - lfo$hyp,
  `Alt1 − Hyp`   = lfo$alt1  - lfo$hyp
) |>
  pivot_longer(-block, names_to = "comparison", values_to = "delta") |>
  mutate(comparison = factor(comparison,
                             levels = c("LogT − Hyp", "LogT − Alt1", "Alt1 − Hyp")))

ann_df <- delta_df |>
  group_by(comparison) |>
  summarise(
    total = sum(delta),
    se    = sd(delta) * sqrt(n()),
    label = sprintf("ΔELPD = %.1f (%.1f SE)", total, abs(total / se)),
    .groups = "drop"
  )

fig3 <- ggplot(delta_df, aes(x = block, y = delta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_col(aes(fill = delta > 0), width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "#5581B4", "FALSE" = "#CC7270")) +
  geom_text(data = ann_df, aes(label = label),
            x = Inf, y = -Inf, hjust = 1.05, vjust = -0.7,
            size = 2.5, family = "Menlo", fontface = "bold", inherit.aes = FALSE) +
  scale_x_continuous(breaks = 6:20) +
  facet_wrap(~ comparison, ncol = 1) +
  labs(x = "Block", 
       y = expression(Delta * "ELPD (LFO-CV)")) +
  theme_classic(base_size = 10, base_family = "Arial") +
  theme(
    strip.background = element_rect(fill = "grey85", color = NA),
    strip.text       = element_text(face = "bold", hjust = 0, size = 9),
    panel.border     = element_rect(color = "grey40", fill = NA, linewidth = 0.5)
  )


ggsave("figure/fig3_lfo_delta.pdf", fig3, 
       width = 5, height = 6.0, units = "in", 
       device = cairo_pdf)


ggsave("figure/fig3_lfo_delta.png", fig3, 
       width = 5, height = 6.0, units = "in", 
       dpi = 600, 
       device = ragg::agg_png)