library(tidyverse)
library(posterior)
library(ggdist)
library(ggtext)
library(grid)

# ============================================================================
# 1. SD for standardization
# ============================================================================

sd_emo  <- sd(df$emots, na.rm = TRUE)
sd_goal <- sd(df$goals, na.rm = TRUE)
sd_perf <- sd(df$perfs, na.rm = TRUE)
sd_disc <- sd(df$perfs - df$goals, na.rm = TRUE)

# ============================================================================
# 2. Parameters actually kept for plotting
# ============================================================================

plot_specs <- tribble(
  ~parameter,        ~label,                              ~panel,                  ~order, ~std_ratio,
  "pop_beta_gp_p",   "&beta;<sup>GP&rarr;P</sup>",        "Effects of the GPD",       1,   sd_disc / sd_perf,
  "pop_beta_gp_g",   "&beta;<sup>GP&rarr;G</sup>",        "Effects of the GPD",       2,   sd_disc / sd_goal,
  "pop_beta_s",      "&beta;<sup>GP&rarr;G|S</sup>",      "Effects of the GPD",       3,   sd_disc / sd_goal,
  "pop_beta_f",      "&beta;<sup>GP&rarr;G|F</sup>",      "Effects of the GPD",       4,   sd_disc / sd_goal,
  "pop_beta_gp_e",   "&beta;<sup>GP&rarr;E</sup>",        "Effects of the GPD",       5,   sd_disc / sd_emo,
  "pop_beta_e_p",    "&beta;<sup>E&rarr;P</sup>",         "Effects of Emotions",      6,   sd_emo  / sd_perf,
  "pop_beta_e_g",    "&beta;<sup>E&rarr;G</sup>",         "Effects of Emotions",      7,   sd_emo  / sd_goal
)

# ============================================================================
# 3. Extract posterior draws
# ============================================================================

get_draws <- function(fit, model_name, params) {
  fit$draws(variables = params, format = "draws_df") |>
    as_tibble() |>
    select(any_of(params)) |>
    pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "value"
    ) |>
    mutate(model = model_name)
}

raw_draws <- bind_rows(
  get_draws(
    fit_hyp, "Hyp",
    c("pop_beta_e_p", "pop_beta_gp_p", "pop_beta_gp_g",
      "pop_beta_e_g", "pop_beta_gp_e")
  ),
  get_draws(
    fit_alt1, "Alt1",
    c("pop_beta_e_p", "pop_beta_gp_p", "pop_beta_s",
      "pop_beta_f", "pop_beta_e_g", "pop_beta_gp_e")
  ),
  get_draws(
    fit_m1, "LogT",
    c("pop_beta_e_p", "pop_beta_gp_p", "pop_beta_s",
      "pop_beta_f", "pop_beta_e_g", "pop_beta_gp_e")
  )
)

# ============================================================================
# 4. Standardize and prepare for plotting
# ============================================================================

draws_std <- raw_draws |>
  inner_join(plot_specs, by = "parameter") |>
  mutate(
    plot_value = value * std_ratio,
    label = factor(
      label,
      levels = plot_specs |>
        arrange(order) |>
        pull(label)
    ),
    panel = factor(
      panel,
      levels = c("Effects of the GPD", "Effects of Emotions")
    ),
    model = factor(model, levels = c("Hyp", "Alt1", "LogT"))
  )

# ============================================================================
# 5. Plot
# ============================================================================

model_fills <- c(
  "Hyp"  = "#E6A817",
  "Alt1" = "#6497B1",
  "LogT" = "#482677"
)

fig4 <- ggplot(draws_std, aes(x = plot_value, y = label, fill = model)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.4
  ) +
  stat_halfeye(
    alpha = 0.75,
    scale = 1,
    normalize = "groups",
    point_size = 0.8,
    linewidth = 0.8,
    .width = 0.95,
    position = position_dodge(width = 0.7),
    color = "black"
  ) +
  scale_fill_manual(
    values = model_fills,
    name = NULL,
    guide = guide_legend(
      direction = "vertical",
      ncol = 1,
      byrow = TRUE
    )
  ) +
  facet_grid(
    rows = vars(panel),
    scales = "free_y",
    space  = "free_y"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 10, base_family = "Arial") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.88, 0.01),
    legend.justification = c(0, 0),
    legend.direction = "vertical",
    
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 8.5, color = "black"),
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(1, "pt"),
    
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    
    axis.ticks.y = element_line(color = "black"),
    axis.text.y = element_markdown(
      hjust = 1,
      size = 10,
      lineheight = 1.05,
      color = "black"
    ),
    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.line = element_blank(),
    
    strip.background = element_rect(fill = "grey85", color = NA),
    strip.text = element_text(
      face = "bold",
      size = 10,
      color = "black",
      margin = margin(b = 5)
    ),
    
    panel.spacing.y = unit(0.6, "lines")
  )

ggsave(
  "figure/fig4_std.pdf",
  fig4,
  width = 6.5,
  height = 7,
  units = "in",
  device = cairo_pdf
)

ggsave(
  "figure/fig4_std.png",
  fig4,
  width = 6.5,
  height = 7,
  units = "in",
  dpi = 300,
  device = ragg::agg_png
)

ggsave(
  "figure/fig4_std.tiff",
  fig4,
  width = 6.5,
  height = 7,
  units = "in",
  dpi = 300,
  device = ragg::agg_tiff,
  compression = "lzw"
)
