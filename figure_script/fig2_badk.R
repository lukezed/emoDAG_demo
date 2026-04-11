library(tidyverse)
library(loo)
library(patchwork)

# ============================================================================
# Panels a–c: plot(loo) returns a ggplot, patch them together
# ============================================================================


make_k_plot <- function(loo_obj, title_label) {
  k <- loo_obj$diagnostics$pareto_k
  df <- tibble(
    obs  = seq_along(k),
    k    = k,
    kcat = case_when(k > 1.0 ~ "k > 1.0", k > 0.7 ~ "k > 0.7", TRUE ~ "OK")
  )
  
  n07 <- sum(k > 0.7); n10 <- sum(k > 1.0); ntot <- length(k)
  ann <- sprintf("k > 0.7: %d (%.1f%%)\nk > 1.0: %d (%.1f%%)",
                 n07, 100*n07/ntot, n10, 100*n10/ntot)
  
  ggplot(df, aes(x = obs, y = k, color = kcat)) +
    geom_point(shape = 3, size = 0.9, stroke = 0.5) +
    geom_hline(yintercept = 0,   linetype = "dotted",  color = "grey60",  linewidth = 0.4) +
    geom_hline(yintercept = 0.7, linetype = "dashed",  color = "#D9A0A0", linewidth = 0.5) +
    geom_hline(yintercept = 1.0, linetype = "solid",   color = "#B04040", linewidth = 0.5) +
    scale_color_manual(
      values = c("OK" = "#6497B1", "k > 0.7" = "#E6A817", "k > 1.0" = "#CC2936"),
      guide  = "none"
    ) +
    annotate("text", x = Inf, y = Inf, label = ann,
             hjust = 1.05, vjust = 1.1, size = 2.5, fontface = "bold", family = "Menlo") +
    labs(x = "Data point", y = "Pareto shape k", title = title_label) +
    theme_classic(base_size = 10, base_family = "Arial") +
    theme(
      plot.title = element_text(size = 9, face = "bold")
    )
}


pa <- make_k_plot(loo_hyp, "Hypothesized")
pb <- make_k_plot(loo_alt1, "Alternate 1")
pc <- make_k_plot(loo_m1, "Log Trial")


fig2 <- (pa | pb) / (pc | plot_spacer()) +
  plot_annotation(
    tag_levels = 'a', 
    tag_prefix = '(', 
    tag_suffix = ')'
  ) +  
  theme(plot.tag = element_text(size = 10, face = "bold", family = "Arial"))


print(fig2)


ggsave("figure/fig2_pareto_k.pdf", fig2, 
       width = 6.5, height = 6.0, units = "in", 
       device = cairo_pdf)

ggsave("figure/fig2_pareto_k.png", fig2, 
       width = 6.5, height = 6.0, units = "in", 
       dpi = 600, 
       device = ragg::agg_png)