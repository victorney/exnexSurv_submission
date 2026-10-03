# Generate only the three figures included in the article.
# Run 03_analyze_article.R first for the simulation summaries.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

export_dir <- "fitted_models/article_exports"
tcga_dir <- "fitted_models/tcga_application"
figure_dir <- "figures"
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

required_files <- c(
  file.path(export_dir, "07_mixed_cohort_operating_characteristics.csv"),
  file.path(export_dir, "04_runtime_trial_level.csv"),
  file.path(tcga_dir, "FigureData_Article_Fig3_TCGA.csv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing input files: ", paste(missing_files, collapse = ", "))
}

read_export <- function(path) {
  read.csv(path, check.names = FALSE)
}

mixed <- read_export(file.path(
  export_dir,
  "07_mixed_cohort_operating_characteristics.csv"
))
runtime_trial <- read_export(file.path(
  export_dir,
  "04_runtime_trial_level.csv"
))
fig_data <- read_export(file.path(
  tcga_dir,
  "FigureData_Article_Fig3_TCGA.csv"
))

method_levels <- c(
  "Complete pooling",
  "No pooling (stratified)",
  "Stan EXNEX (NUTS)",
  "exnexSurv (Gibbs DA)"
)
method_colors <- c(
  "Complete pooling" = "#D55E00",
  "No pooling (stratified)" = "#E69F00",
  "Stan EXNEX (NUTS)" = "#009E73",
  "exnexSurv (Gibbs DA)" = "#0072B2"
)
scenario_levels <- c(
  "Homogeneous",
  "Mixed efficacy",
  "Complete heterogeneity"
)

paper_theme <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.6),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey25", hjust = 0),
      plot.margin = margin(8, 10, 16, 8)
    )
}

wrap_labels <- function(x, width = 12) {
  vapply(
    x,
    function(s) paste(strwrap(s, width = width), collapse = "\n"),
    character(1)
  )
}

# Main Figure 2: mixed-efficacy cohort performance.
mixed_primary <- mixed %>%
  filter(scenario == "mixed", population == "qualified_trials") %>%
  mutate(
    cohort = factor(
      cohort,
      levels = c(
        "Responsive (baskets 1-6)",
        "Resistant outliers (baskets 7-9)"
      )
    ),
    method = factor(method, levels = method_levels)
  )

p2_bias <- ggplot(mixed_primary, aes(cohort, bias, fill = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  geom_col(
    position = position_dodge(0.78),
    width = 0.7,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = method_colors, name = "Method", drop = FALSE) +
  scale_x_discrete(labels = wrap_labels) +
  labs(title = "(A) Bias", x = NULL, y = "Bias on log-time scale") +
  paper_theme() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

p2_rmse <- ggplot(mixed_primary, aes(cohort, rmse, fill = method)) +
  geom_col(
    position = position_dodge(0.78),
    width = 0.7,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = method_colors, name = "Method", drop = FALSE) +
  scale_x_discrete(labels = wrap_labels) +
  labs(title = "(B) RMSE", x = NULL, y = "RMSE on log-time scale") +
  paper_theme() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

p2_coverage <- ggplot(mixed_primary, aes(cohort, coverage_95, fill = method)) +
  geom_col(
    position = position_dodge(0.78),
    width = 0.7,
    color = "black",
    linewidth = 0.2
  ) +
  geom_errorbar(
    aes(ymin = coverage_95_lower_mc, ymax = coverage_95_upper_mc),
    position = position_dodge(0.78),
    width = 0.18,
    linewidth = 0.35
  ) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "grey30") +
  scale_fill_manual(values = method_colors, name = "Method", drop = FALSE) +
  scale_x_discrete(labels = wrap_labels) +
  scale_y_continuous(limits = c(0, 105)) +
  labs(title = "(C) 95% interval coverage", x = NULL, y = "Coverage (%)") +
  paper_theme() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

fig2 <- (p2_bias + p2_rmse + p2_coverage) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")
ggsave(
  file.path(figure_dir, "Main_Fig2_Mixed_Cohorts.pdf"),
  fig2,
  width = 13,
  height = 5.2,
  units = "in",
  device = "pdf",
  useDingbats = FALSE
)

# Main Figure 3: paired computational performance.
comparison_methods <- c("exnexSurv (Gibbs DA)", "Stan EXNEX (NUTS)")
comparison_runtime <- runtime_trial %>%
  filter(
    method %in% comparison_methods,
    scenario %in% c("homogeneous", "mixed", "heterogeneous")
  ) %>%
  mutate(
    method = factor(method, levels = comparison_methods),
    scenario_label = factor(scenario_label, levels = scenario_levels)
  )

runtime_plot_data <- comparison_runtime %>%
  select(scenario, scenario_label, rep_id, method, runtime_seconds)
ess_plot_data <- comparison_runtime %>%
  select(scenario, scenario_label, rep_id, method, ess_per_second)
speedup_plot_data <- comparison_runtime %>%
  select(scenario, scenario_label, rep_id, method, runtime_seconds) %>%
  pivot_wider(names_from = method, values_from = runtime_seconds) %>%
  transmute(
    scenario,
    scenario_label,
    rep_id,
    speedup_stan_over_gibbs = `Stan EXNEX (NUTS)` / `exnexSurv (Gibbs DA)`
  )

p3_runtime <- ggplot(
  runtime_plot_data,
  aes(scenario_label, runtime_seconds, fill = method)
) +
  geom_boxplot(outlier.size = 0.25, outlier.alpha = 0.15, width = 0.7) +
  scale_y_log10() +
  scale_fill_manual(values = method_colors, name = "Method", drop = FALSE) +
  scale_x_discrete(labels = wrap_labels) +
  labs(title = "(A) Runtime per trial", x = NULL, y = "Seconds (log scale)") +
  paper_theme() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

p3_ess <- ggplot(
  ess_plot_data,
  aes(scenario_label, ess_per_second, fill = method)
) +
  geom_boxplot(outlier.size = 0.25, outlier.alpha = 0.15, width = 0.7) +
  scale_y_log10() +
  scale_fill_manual(values = method_colors, name = "Method", drop = FALSE) +
  scale_x_discrete(labels = wrap_labels) +
  labs(
    title = "(B) Effective sample size per second",
    x = NULL,
    y = "ESS/sec (log scale)"
  ) +
  paper_theme() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

p3_speedup <- ggplot(
  speedup_plot_data,
  aes(scenario_label, speedup_stan_over_gibbs)
) +
  geom_boxplot(
    outlier.size = 0.25,
    outlier.alpha = 0.15,
    width = 0.6,
    fill = "grey75"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey30") +
  scale_x_discrete(labels = wrap_labels) +
  labs(title = "(C) Stan/Gibbs runtime ratio", x = NULL, y = "Speed-up ratio") +
  paper_theme() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none"
  )

fig3 <- (p3_runtime + p3_ess + p3_speedup) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")
ggsave(
  file.path(figure_dir, "Main_Fig3_Computational.pdf"),
  fig3,
  width = 13,
  height = 5.2,
  units = "in",
  device = "pdf",
  useDingbats = FALSE
)

# Article Figure 3: TCGA application.
paper_theme_tcga <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.6),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey25", hjust = 0),
      plot.margin = margin(6, 8, 6, 6)
    )
}
method_colors_tcga <- c(
  "exnexSurv (Gibbs DA)" = "#0072B2",
  "Stan EXNEX (NUTS)" = "#009E73",
  "No pooling (stratified)" = "#E69F00"
)

intercepts <- fig_data[fig_data$panel == "A_baskets", ]
pooled <- fig_data[fig_data$panel == "A_pooled_band", ]
basket_order <- intercepts %>%
  filter(method == "exnexSurv (Gibbs DA)") %>%
  arrange(estimate) %>%
  pull(basket)
intercepts$basket <- factor(intercepts$basket, levels = basket_order)

panel_a <- ggplot() +
  geom_rect(
    data = pooled,
    aes(xmin = 0.4, xmax = 9.6, ymin = q025, ymax = q975),
    fill = "grey70",
    alpha = 0.25
  ) +
  geom_hline(
    data = pooled,
    aes(yintercept = estimate),
    linetype = "dashed",
    color = "grey30"
  ) +
  geom_pointrange(
    data = intercepts,
    aes(basket, estimate, ymin = q025, ymax = q975, color = method),
    position = position_dodge(width = 0.7),
    size = 0.35,
    linewidth = 0.5
  ) +
  scale_color_manual(
    values = method_colors_tcga,
    name = "Method",
    drop = FALSE
  ) +
  labs(
    title = "(A) Basket-specific log-survival intercepts",
    subtitle = "Gray band: Complete-pooling estimate (95% CrI)",
    x = NULL,
    y = "Posterior mean (95% CrI)"
  ) +
  paper_theme_tcga()

gibbs <- fig_data %>%
  filter(panel == "B_shrinkage", method == "exnexSurv (Gibbs DA)") %>%
  select(basket, estimate, q025, q975, n_events = n)
no_pool <- fig_data %>%
  filter(panel == "A_baskets", method == "No pooling (stratified)") %>%
  select(basket, np_est = estimate, np_q025 = q025, np_q975 = q975)
shrinkage <- left_join(gibbs, no_pool, by = "basket")

panel_b <- ggplot(shrinkage, aes(np_est, estimate)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey30") +
  geom_errorbar(
    aes(y = estimate, xmin = np_q025, xmax = np_q975),
    orientation = "y",
    width = 0,
    color = "grey40"
  ) +
  geom_errorbar(
    aes(ymin = q025, ymax = q975),
    width = 0,
    color = "grey40"
  ) +
  geom_point(aes(size = n_events), color = "#0072B2") +
  geom_text(aes(label = basket), vjust = -1.1, size = 2.6) +
  scale_size_continuous(range = c(1.5, 5)) +
  labs(
    title = "(B) Shrinkage of basket effects under EXNEX",
    x = "No pooling posterior mean",
    y = "exnexSurv posterior mean",
    size = "Events"
  ) +
  paper_theme_tcga()

tcga_figure <- panel_a +
  panel_b +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(
  file.path(figure_dir, "Article_Fig3_TCGA.pdf"),
  tcga_figure,
  width = 11,
  height = 4.5,
  units = "in",
  device = "pdf",
  useDingbats = FALSE
)

cat("Wrote the three article figures to ", figure_dir, ".\n", sep = "")
