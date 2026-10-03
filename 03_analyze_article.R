# Derive the article evidence exports from the archived fit summaries.
#
# Reads fitted_models/scenarios/results_scenario_*.rds and the TCGA application
# summaries, and writes the compact CSVs in fitted_models/article_exports. It
# does not fit any model.
#
# Usage: Rscript 03_analyze_article.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

result_dir <- "fitted_models/scenarios"
export_dir <- "fitted_models/article_exports"
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

scenario_files <- list.files(
  result_dir,
  pattern = "^results_scenario_.*\\.rds$",
  full.names = TRUE
)
# The divergent-position checkpoints are analyzed by 08_divergent_position.R.
scenario_files <- scenario_files[
  !grepl("mixed_divlarge", basename(scenario_files))
]
if (length(scenario_files) == 0L) {
  stop("No archived scenario results found in ", result_dir, ".")
}

df_raw <- bind_rows(lapply(scenario_files, readRDS))
if (!"target_cens" %in% names(df_raw)) {
  df_raw$target_cens <- NA_real_
}
if (!"observed_cens" %in% names(df_raw)) {
  df_raw$observed_cens <- NA_real_
}

required_columns <- c(
  "scenario",
  "rep_id",
  "group",
  "theta_true",
  "n_total",
  "n_events",
  "time_exnex",
  "time_stan_exnex",
  "time_stan_pooled",
  "time_stan_unp",
  "exnex_mean",
  "exnex_sd",
  "exnex_q025",
  "exnex_q975",
  "exnex_rhat",
  "exnex_ess",
  "stan_exnex.mean",
  "stan_exnex.sd",
  "stan_exnex.q025",
  "stan_exnex.q975",
  "stan_exnex.rhat",
  "stan_exnex.ess",
  "stan_pooled.mean",
  "stan_pooled.sd",
  "stan_pooled.q025",
  "stan_pooled.q975",
  "stan_pooled.rhat",
  "stan_pooled.ess",
  "stan_unpooled.mean",
  "stan_unpooled.sd",
  "stan_unpooled.q025",
  "stan_unpooled.q975",
  "stan_unpooled.rhat",
  "stan_unpooled.ess"
)
missing_columns <- setdiff(required_columns, names(df_raw))
if (length(missing_columns) > 0L) {
  stop(
    "Archived results are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

df_raw <- df_raw %>%
  mutate(
    scenario = as.character(scenario),
    rep_id = as.integer(rep_id),
    group = as.integer(group),
    n_total = as.numeric(n_total),
    n_events = as.numeric(n_events)
  )

scenario_order <- c(
  "homogeneous",
  "mixed",
  "heterogeneous",
  "weibull_misspec",
  "hetvar_misspec",
  "t_misspec",
  "infcens_misspec",
  "tcga_calib"
)
scenario_labels <- c(
  homogeneous = "Homogeneous",
  mixed = "Mixed efficacy",
  heterogeneous = "Complete heterogeneity",
  weibull_misspec = "Gumbel (extreme-value) misspecification",
  hetvar_misspec = "Heteroskedastic variances",
  t_misspec = "Student-t errors",
  infcens_misspec = "Informative censoring",
  tcga_calib = "TCGA-calibrated"
)
scenario_definitions <- c(
  homogeneous = "Nine baskets generated around a common positive mean with low between-basket heterogeneity.",
  mixed = "Six responsive baskets around a positive mean and three resistant outlier baskets.",
  heterogeneous = "Nine fixed basket effects spanning a broad heterogeneous range.",
  weibull_misspec = "Mixed-efficacy effects with Gumbel (extreme-value) errors on the log-time scale.",
  hetvar_misspec = "Mixed-efficacy effects with basket-specific residual variances.",
  t_misspec = "Mixed-efficacy effects with scaled Student-t errors (5 degrees of freedom).",
  infcens_misspec = "Mixed-efficacy effects with informative censoring.",
  tcga_calib = "TCGA-calibrated fixed effects, basket sizes, and censoring targets."
)

unknown_scenarios <- setdiff(unique(df_raw$scenario), scenario_order)
if (length(unknown_scenarios) > 0L) {
  stop("Unrecognized scenarios: ", paste(unknown_scenarios, collapse = ", "))
}

df_raw <- df_raw %>%
  mutate(
    scenario_label = unname(scenario_labels[scenario]),
    cohort = case_when(
      scenario %in%
        c(
          "mixed",
          "weibull_misspec",
          "hetvar_misspec",
          "t_misspec",
          "infcens_misspec"
        ) &
        group <= 6L ~ "Responsive (baskets 1-6)",
      scenario %in%
        c(
          "mixed",
          "weibull_misspec",
          "hetvar_misspec",
          "t_misspec",
          "infcens_misspec"
        ) &
        group >= 7L ~ "Resistant outliers (baskets 7-9)",
      TRUE ~ NA_character_
    )
  )

trial_censoring <- df_raw %>%
  group_by(scenario, scenario_label, rep_id) %>%
  summarise(
    target_cens = coalesce(first(target_cens), 0.30),
    n_total_trial = sum(n_total),
    n_events_trial = sum(n_events),
    observed_cens = 1 - n_events_trial / n_total_trial,
    .groups = "drop"
  )
if (any(!is.finite(trial_censoring$observed_cens))) {
  stop("Archived results contain non-finite censoring proportions.")
}

df_raw <- df_raw %>%
  select(-any_of(c("target_cens", "observed_cens"))) %>%
  left_join(
    trial_censoring %>% select(scenario, rep_id, target_cens, observed_cens),
    by = c("scenario", "rep_id")
  )

method_specs <- tibble(
  method_id = c("gibbs_exnex", "stan_exnex", "stan_pooled", "stan_unpooled"),
  method = c(
    "exnexSurv (Gibbs DA)",
    "Stan EXNEX (NUTS)",
    "Complete pooling",
    "No pooling (stratified)"
  ),
  estimate = c(
    "exnex_mean",
    "stan_exnex.mean",
    "stan_pooled.mean",
    "stan_unpooled.mean"
  ),
  posterior_sd = c(
    "exnex_sd",
    "stan_exnex.sd",
    "stan_pooled.sd",
    "stan_unpooled.sd"
  ),
  q025 = c(
    "exnex_q025",
    "stan_exnex.q025",
    "stan_pooled.q025",
    "stan_unpooled.q025"
  ),
  q975 = c(
    "exnex_q975",
    "stan_exnex.q975",
    "stan_pooled.q975",
    "stan_unpooled.q975"
  ),
  rhat = c(
    "exnex_rhat",
    "stan_exnex.rhat",
    "stan_pooled.rhat",
    "stan_unpooled.rhat"
  ),
  ess = c(
    "exnex_ess",
    "stan_exnex.ess",
    "stan_pooled.ess",
    "stan_unpooled.ess"
  ),
  runtime = c(
    "time_exnex",
    "time_stan_exnex",
    "time_stan_pooled",
    "time_stan_unp"
  )
)

euler_gamma <- -digamma(1)
sigma_w <- 1.2 * sqrt(6 / pi^2)

long_results <- bind_rows(lapply(seq_len(nrow(method_specs)), function(i) {
  spec <- method_specs[i, ]
  df_raw %>%
    transmute(
      scenario,
      scenario_label,
      rep_id,
      group,
      cohort,
      target_cens,
      observed_cens,
      n_total,
      n_events,
      theta_true,
      estimate = .data[[spec$estimate]],
      posterior_sd = .data[[spec$posterior_sd]],
      q025 = .data[[spec$q025]],
      q975 = .data[[spec$q975]],
      rhat = .data[[spec$rhat]],
      bulk_ess = .data[[spec$ess]],
      runtime_seconds = .data[[spec$runtime]],
      method_id = spec$method_id,
      method = spec$method
    ) %>%
    mutate(
      theta_estimand = ifelse(
        scenario == "weibull_misspec",
        theta_true + euler_gamma * sigma_w,
        theta_true
      ),
      error = estimate - theta_estimand,
      interval_width = q975 - q025,
      covered_95 = theta_estimand >= q025 & theta_estimand <= q975,
      finite_summary = is.finite(estimate) & is.finite(q025) & is.finite(q975),
      ordered_interval = is.finite(q025) & is.finite(q975) & q025 <= q975,
      rhat_lt_1_01 = is.finite(rhat) & rhat < 1.01,
      rhat_lt_1_05 = is.finite(rhat) & rhat < 1.05,
      ess_ge_100 = is.finite(bulk_ess) & bulk_ess >= 100,
      ess_ge_400 = is.finite(bulk_ess) & bulk_ess >= 400,
      qualified = finite_summary & ordered_interval & rhat_lt_1_01 & ess_ge_400
    )
}))

trial_status <- long_results %>%
  group_by(scenario, scenario_label, rep_id, method_id, method) %>%
  summarise(
    n_baskets = n_distinct(group),
    qualified_trial = n_baskets == 9L & all(qualified),
    .groups = "drop"
  )
long_results <- long_results %>%
  left_join(
    trial_status %>% select(scenario, rep_id, method_id, qualified_trial),
    by = c("scenario", "rep_id", "method_id")
  )

long_all <- long_results %>% mutate(population = "all_fits")
long_qualified <- long_results %>%
  filter(qualified_trial) %>%
  mutate(population = "qualified_trials")
joint_trials <- trial_status %>%
  filter(method_id %in% c("gibbs_exnex", "stan_exnex")) %>%
  group_by(scenario, rep_id) %>%
  summarise(
    joint_qualified = n() == 2L & all(qualified_trial),
    .groups = "drop"
  )
long_joint <- long_results %>%
  filter(method_id %in% c("gibbs_exnex", "stan_exnex")) %>%
  left_join(joint_trials, by = c("scenario", "rep_id")) %>%
  filter(joint_qualified) %>%
  mutate(population = "jointly_qualified")
long_populations <- bind_rows(long_all, long_qualified, long_joint) %>%
  filter(is.finite(estimate))

safe_sd <- function(x) {
  if (sum(is.finite(x)) > 1L) sd(x, na.rm = TRUE) else NA_real_
}
safe_quantile <- function(x, p) {
  if (sum(is.finite(x)) > 0L) {
    as.numeric(quantile(x, p, na.rm = TRUE, names = FALSE))
  } else {
    NA_real_
  }
}
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) > 1L) cor(x[ok], y[ok]) else NA_real_
}
safe_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) > 1L && safe_sd(x[ok]) > 0) {
    as.numeric(coef(lm(y[ok] ~ x[ok]))[[2L]])
  } else {
    NA_real_
  }
}

summarize_operating <- function(data, by_group = FALSE) {
  summary_vars <- c(
    "scenario",
    "scenario_label",
    "method_id",
    "method",
    "population",
    "cohort",
    if (by_group) "group"
  )
  trial_vars <- c(
    "scenario",
    "scenario_label",
    "rep_id",
    "method_id",
    "method",
    "population",
    "cohort",
    if (by_group) "group"
  )

  trial_metrics <- data %>%
    group_by(across(all_of(trial_vars))) %>%
    summarise(
      trial_bias = mean(error, na.rm = TRUE),
      trial_mse = mean(error^2, na.rm = TRUE),
      trial_coverage = mean(covered_95, na.rm = TRUE),
      trial_width = mean(interval_width, na.rm = TRUE),
      .groups = "drop"
    )
  trial_summary <- trial_metrics %>%
    group_by(across(all_of(summary_vars))) %>%
    summarise(
      n_trials = n_distinct(rep_id),
      mcse_bias = safe_sd(trial_bias) / sqrt(n()),
      mean_trial_mse = mean(trial_mse, na.rm = TRUE),
      mcse_trial_mse = safe_sd(trial_mse) / sqrt(n()),
      mcse_coverage = safe_sd(trial_coverage) / sqrt(n()) * 100,
      mcse_interval_width = safe_sd(trial_width) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      rmse = sqrt(mean_trial_mse),
      mcse_rmse = ifelse(rmse > 0, mcse_trial_mse / (2 * rmse), NA_real_)
    )

  data %>%
    group_by(across(all_of(summary_vars))) %>%
    summarise(
      n_basket_evaluations = n(),
      n_qualified_trials = ifelse(
        first(population) == "qualified_trials",
        n_distinct(rep_id),
        NA_integer_
      ),
      bias = mean(error, na.rm = TRUE),
      empirical_sd = safe_sd(estimate),
      error_sd = safe_sd(error),
      coverage_95 = mean(covered_95, na.rm = TRUE) * 100,
      mean_interval_width = mean(interval_width, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(trial_summary, by = summary_vars) %>%
    mutate(
      coverage_95_lower_mc = pmax(0, coverage_95 - 1.96 * mcse_coverage),
      coverage_95_upper_mc = pmin(100, coverage_95 + 1.96 * mcse_coverage)
    )
}

global_operating <- summarize_operating(
  long_populations %>% mutate(cohort = NA_character_)
)
cohort_operating <- summarize_operating(
  long_populations %>%
    filter(
      scenario %in%
        c(
          "mixed",
          "weibull_misspec",
          "hetvar_misspec",
          "t_misspec",
          "infcens_misspec"
        ),
      !is.na(cohort)
    )
)

trial_censoring_summary <- trial_censoring %>%
  group_by(scenario, scenario_label) %>%
  summarise(
    n_trials = n(),
    target_censoring = mean(target_cens),
    mean_observed_censoring = mean(observed_cens),
    sd_observed_censoring = safe_sd(observed_cens),
    median_observed_censoring = median(observed_cens),
    observed_censoring_p025 = safe_quantile(observed_cens, 0.025),
    observed_censoring_p975 = safe_quantile(observed_cens, 0.975),
    mean_patients_per_trial = mean(n_total_trial),
    mean_events_per_trial = mean(n_events_trial),
    .groups = "drop"
  )
study_design <- trial_censoring_summary %>%
  mutate(
    n_baskets = 9L,
    total_patients_per_trial = round(mean_patients_per_trial),
    scenario_definition = unname(scenario_definitions[scenario])
  ) %>%
  select(
    scenario,
    scenario_label,
    scenario_definition,
    n_trials,
    n_baskets,
    total_patients_per_trial,
    target_censoring,
    mean_observed_censoring,
    sd_observed_censoring,
    mean_patients_per_trial,
    mean_events_per_trial
  )

runtime_trial <- df_raw %>%
  group_by(scenario, scenario_label, rep_id) %>%
  summarise(
    target_cens = first(target_cens),
    observed_cens = first(observed_cens),
    time_gibbs = first(time_exnex),
    time_stan_exnex = first(time_stan_exnex),
    time_stan_pooled = first(time_stan_pooled),
    time_stan_unpooled = first(time_stan_unp),
    gibbs_mean_ess = mean(exnex_ess, na.rm = TRUE),
    stan_exnex_mean_ess = mean(stan_exnex.ess, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    speedup_stan_over_gibbs = time_stan_exnex / time_gibbs,
    gibbs_ess_per_second = gibbs_mean_ess / time_gibbs,
    stan_exnex_ess_per_second = stan_exnex_mean_ess / time_stan_exnex
  )
runtime_long <- bind_rows(
  runtime_trial %>%
    transmute(
      scenario,
      scenario_label,
      rep_id,
      method = "exnexSurv (Gibbs DA)",
      runtime_seconds = time_gibbs,
      mean_ess = gibbs_mean_ess,
      ess_per_second = gibbs_ess_per_second
    ),
  runtime_trial %>%
    transmute(
      scenario,
      scenario_label,
      rep_id,
      method = "Stan EXNEX (NUTS)",
      runtime_seconds = time_stan_exnex,
      mean_ess = stan_exnex_mean_ess,
      ess_per_second = stan_exnex_ess_per_second
    ),
  runtime_trial %>%
    transmute(
      scenario,
      scenario_label,
      rep_id,
      method = "Complete pooling",
      runtime_seconds = time_stan_pooled,
      mean_ess = NA_real_,
      ess_per_second = NA_real_
    ),
  runtime_trial %>%
    transmute(
      scenario,
      scenario_label,
      rep_id,
      method = "No pooling (stratified)",
      runtime_seconds = time_stan_unpooled,
      mean_ess = NA_real_,
      ess_per_second = NA_real_
    )
)
runtime_summary <- runtime_long %>%
  group_by(scenario, scenario_label, method) %>%
  summarise(
    n_trials = sum(!is.na(runtime_seconds)),
    runtime_mean_seconds = mean(runtime_seconds, na.rm = TRUE),
    runtime_sd_seconds = safe_sd(runtime_seconds),
    runtime_median_seconds = median(runtime_seconds, na.rm = TRUE),
    runtime_iqr_seconds = IQR(runtime_seconds, na.rm = TRUE),
    runtime_p025_seconds = safe_quantile(runtime_seconds, 0.025),
    runtime_p975_seconds = safe_quantile(runtime_seconds, 0.975),
    mean_ess = mean(mean_ess, na.rm = TRUE),
    mean_ess_per_second = mean(ess_per_second, na.rm = TRUE),
    median_ess_per_second = median(ess_per_second, na.rm = TRUE),
    .groups = "drop"
  )
speedup_summary <- runtime_trial %>%
  group_by(scenario, scenario_label) %>%
  summarise(
    n_trials = sum(!is.na(speedup_stan_over_gibbs)),
    median_speedup_stan_over_gibbs = median(
      speedup_stan_over_gibbs,
      na.rm = TRUE
    ),
    mean_speedup_stan_over_gibbs = mean(speedup_stan_over_gibbs, na.rm = TRUE),
    speedup_p025 = safe_quantile(speedup_stan_over_gibbs, 0.025),
    speedup_p975 = safe_quantile(speedup_stan_over_gibbs, 0.975),
    .groups = "drop"
  )

agreement <- df_raw %>%
  transmute(
    scenario,
    scenario_label,
    rep_id,
    group,
    gibbs_mean = exnex_mean,
    stan_mean = stan_exnex.mean
  ) %>%
  left_join(
    long_results %>%
      filter(method_id == "gibbs_exnex") %>%
      select(scenario, rep_id, group, gibbs_qualified = qualified),
    by = c("scenario", "rep_id", "group")
  ) %>%
  left_join(
    long_results %>%
      filter(method_id == "stan_exnex") %>%
      select(scenario, rep_id, group, stan_qualified = qualified),
    by = c("scenario", "rep_id", "group")
  ) %>%
  mutate(
    joint_qualified = gibbs_qualified & stan_qualified,
    difference = gibbs_mean - stan_mean,
    absolute_difference = abs(difference),
    average = (gibbs_mean + stan_mean) / 2
  )
agreement_summary <- function(data, comparison_level, population) {
  data %>%
    filter(is.finite(absolute_difference)) %>%
    group_by(scenario, scenario_label) %>%
    summarise(
      comparison_level = comparison_level,
      population = population,
      n_basket_comparisons = n(),
      n_trials = n_distinct(rep_id),
      mean_absolute_difference = mean(absolute_difference),
      median_absolute_difference = median(absolute_difference),
      p95_absolute_difference = safe_quantile(absolute_difference, 0.95),
      p99_absolute_difference = safe_quantile(absolute_difference, 0.99),
      max_absolute_difference = max(absolute_difference),
      proportion_absolute_difference_gt_0_05 = mean(
        absolute_difference > 0.05
      ) *
        100,
      proportion_absolute_difference_gt_0_10 = mean(
        absolute_difference > 0.10
      ) *
        100,
      proportion_absolute_difference_gt_1 = mean(absolute_difference > 1) * 100,
      pearson_correlation = safe_cor(stan_mean, gibbs_mean),
      identity_regression_slope = safe_slope(stan_mean, gibbs_mean),
      .groups = "drop"
    )
}
agreement_trial <- agreement %>%
  group_by(scenario, scenario_label, rep_id) %>%
  summarise(
    joint_qualified = all(joint_qualified),
    gibbs_mean = mean(gibbs_mean, na.rm = TRUE),
    stan_mean = mean(stan_mean, na.rm = TRUE),
    absolute_difference = abs(gibbs_mean - stan_mean),
    .groups = "drop"
  )
agreement_summary_table <- bind_rows(
  agreement_summary(agreement, "basket", "all_finite_fits"),
  agreement_summary(
    filter(agreement, joint_qualified),
    "basket",
    "jointly_qualified"
  ),
  agreement_summary(agreement_trial, "trial_mean", "all_finite_fits"),
  agreement_summary(
    filter(agreement_trial, joint_qualified),
    "trial_mean",
    "jointly_qualified"
  )
)

paired_data <- agreement %>%
  filter(is.finite(gibbs_mean), is.finite(stan_mean))

paired_basket <- paired_data %>%
  group_by(scenario, scenario_label) %>%
  summarise(
    level = "basket",
    n = n(),
    both = sum(gibbs_qualified & stan_qualified),
    gibbs_only = sum(gibbs_qualified & !stan_qualified),
    stan_only = sum(!gibbs_qualified & stan_qualified),
    neither = sum(!gibbs_qualified & !stan_qualified),
    .groups = "drop"
  )
trial_quals <- paired_data %>%
  group_by(scenario, scenario_label, rep_id) %>%
  summarise(
    gibbs_trial = all(gibbs_qualified),
    stan_trial = all(stan_qualified),
    .groups = "drop"
  )
paired_trial <- trial_quals %>%
  group_by(scenario, scenario_label) %>%
  summarise(
    level = "trial",
    n = n(),
    both = sum(gibbs_trial & stan_trial),
    gibbs_only = sum(gibbs_trial & !stan_trial),
    stan_only = sum(!gibbs_trial & stan_trial),
    neither = sum(!gibbs_trial & !stan_trial),
    .groups = "drop"
  )
paired <- bind_rows(paired_basket, paired_trial) %>%
  mutate(
    both_pct = 100 * both / n,
    gibbs_only_pct = 100 * gibbs_only / n,
    stan_only_pct = 100 * stan_only / n,
    neither_pct = 100 * neither / n
  )

variance_groups <- tibble(
  group = 1:9,
  variance_group = case_when(
    group <= 3L ~ "Low",
    group <= 6L ~ "Medium",
    TRUE ~ "High"
  )
)
hetvar_groups <- long_populations %>%
  filter(scenario == "hetvar_misspec") %>%
  left_join(variance_groups, by = "group") %>%
  group_by(scenario_label, method, population, variance_group) %>%
  summarise(
    n_basket_evaluations = n(),
    bias = mean(error, na.rm = TRUE),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    coverage_95 = mean(covered_95, na.rm = TRUE) * 100,
    mean_interval_width = mean(interval_width, na.rm = TRUE),
    .groups = "drop"
  )

extreme_discrepancies <- long_results %>%
  filter(method_id %in% c("gibbs_exnex", "stan_exnex")) %>%
  select(
    scenario,
    scenario_label,
    rep_id,
    group,
    n_total,
    n_events,
    method_id,
    estimate,
    rhat,
    bulk_ess,
    qualified
  ) %>%
  pivot_wider(
    id_cols = c(scenario, scenario_label, rep_id, group, n_total, n_events),
    names_from = method_id,
    values_from = c(estimate, rhat, bulk_ess, qualified)
  ) %>%
  mutate(
    absolute_difference = abs(estimate_gibbs_exnex - estimate_stan_exnex),
    observed_cens = 1 - n_events / n_total,
    extreme = absolute_difference > 1,
    joint_qualified = qualified_gibbs_exnex & qualified_stan_exnex
  ) %>%
  filter(extreme) %>%
  group_by(scenario) %>%
  arrange(desc(absolute_difference), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup() %>%
  select(
    scenario,
    scenario_label,
    rep_id,
    group,
    n_total,
    n_events,
    observed_cens,
    gibbs_mean = estimate_gibbs_exnex,
    stan_mean = estimate_stan_exnex,
    absolute_difference,
    gibbs_rhat = rhat_gibbs_exnex,
    stan_rhat = rhat_stan_exnex,
    gibbs_ess = bulk_ess_gibbs_exnex,
    stan_ess = bulk_ess_stan_exnex,
    gibbs_qualified = qualified_gibbs_exnex,
    stan_qualified = qualified_stan_exnex,
    joint_qualified
  )

ess_index <- runtime_long %>%
  filter(method %in% c("exnexSurv (Gibbs DA)", "Stan EXNEX (NUTS)")) %>%
  mutate(
    method_id = ifelse(
      method == "exnexSurv (Gibbs DA)",
      "gibbs_exnex",
      "stan_exnex"
    )
  ) %>%
  left_join(
    long_results %>%
      filter(method_id %in% c("gibbs_exnex", "stan_exnex")) %>%
      group_by(scenario, scenario_label, rep_id, method_id) %>%
      summarise(
        min_basket_ess = if (any(is.finite(bulk_ess))) {
          min(bulk_ess[is.finite(bulk_ess)])
        } else {
          NA_real_
        },
        qualified_trial = first(qualified_trial),
        .groups = "drop"
      ),
    by = c("scenario", "scenario_label", "rep_id", "method_id")
  ) %>%
  mutate(ess_normalized_seconds = runtime_seconds * 400 / min_basket_ess) %>%
  group_by(scenario, rep_id) %>%
  mutate(jointly_qualified = all(qualified_trial == TRUE)) %>%
  ungroup()
ess_index_summary <- ess_index %>%
  filter(jointly_qualified, is.finite(ess_normalized_seconds)) %>%
  group_by(scenario, scenario_label, method, method_id) %>%
  summarise(
    n_trials = n_distinct(rep_id),
    median_ess_normalized_seconds = median(ess_normalized_seconds),
    p25_ess_normalized_seconds = quantile(
      ess_normalized_seconds,
      0.25,
      names = FALSE
    ),
    p75_ess_normalized_seconds = quantile(
      ess_normalized_seconds,
      0.75,
      names = FALSE
    ),
    .groups = "drop"
  )

# Accuracy-matched index for the TCGA application itself (Table 3). The two
# CSVs are written by 05_tcga_application.R.
app_runtime <- read.csv("fitted_models/tcga_application/04_runtime.csv")
app_diagnostics <- read.csv(
  "fitted_models/tcga_application/03_mcmc_diagnostics.csv"
)
app_index <- app_diagnostics %>%
  filter(
    parameter == "theta",
    method_id %in% c("gibbs_exnex", "stan_exnex")
  ) %>%
  group_by(method_id, method) %>%
  summarise(min_basket_ess = min(bulk_ess), .groups = "drop") %>%
  left_join(
    app_runtime %>% select(method_id, runtime_seconds = seconds),
    by = "method_id"
  ) %>%
  mutate(
    scenario = "tcga_application",
    scenario_label = "TCGA application",
    n_trials = 1L,
    ess_normalized_seconds = runtime_seconds * 400 / min_basket_ess
  ) %>%
  transmute(
    scenario,
    scenario_label,
    method,
    method_id,
    n_trials,
    median_ess_normalized_seconds = ess_normalized_seconds,
    p25_ess_normalized_seconds = ess_normalized_seconds,
    p75_ess_normalized_seconds = ess_normalized_seconds
  )
ess_index_summary <- bind_rows(ess_index_summary, app_index)

# TCGA-calibrated qualification: the Gibbs sampler qualifies on 758 of the
# 1,000 fitted replicates; on the 100 trials where the Stan comparator was
# run, Gibbs qualifies on 71, Stan on 99, and both on 70.
calib <- readRDS("fitted_models/scenarios/results_scenario_tcga_calib.rds")
calib_trials <- calib %>%
  mutate(
    gibbs_basket = is.finite(exnex_mean) &
      is.finite(exnex_q025) &
      is.finite(exnex_q975) &
      exnex_q025 <= exnex_q975 &
      is.finite(exnex_rhat) &
      exnex_rhat < 1.01 &
      is.finite(exnex_ess) &
      exnex_ess >= 400,
    stan_basket = is.finite(stan_exnex.mean) &
      is.finite(stan_exnex.q025) &
      is.finite(stan_exnex.q975) &
      stan_exnex.q025 <= stan_exnex.q975 &
      is.finite(stan_exnex.rhat) &
      stan_exnex.rhat < 1.01 &
      is.finite(stan_exnex.ess) &
      stan_exnex.ess >= 400
  ) %>%
  group_by(rep_id) %>%
  summarise(
    has_stan = any(is.finite(stan_exnex.mean)),
    gibbs_qualified = sum(gibbs_basket) == n(),
    stan_qualified = any(is.finite(stan_exnex.mean)) & sum(stan_basket) == n(),
    .groups = "drop"
  )
calib_qualification <- bind_rows(
  data.frame(
    scenario = "tcga_calib",
    population = "all_replicates",
    method = "gibbs_exnex",
    n_trials = nrow(calib_trials),
    n_qualified = sum(calib_trials$gibbs_qualified)
  ),
  data.frame(
    scenario = "tcga_calib",
    population = "stan_subset",
    method = c("gibbs_exnex", "stan_exnex", "jointly_qualified"),
    n_trials = 100L,
    n_qualified = c(
      sum(calib_trials$gibbs_qualified[calib_trials$rep_id <= 100]),
      sum(calib_trials$stan_qualified[calib_trials$rep_id <= 100]),
      sum(
        calib_trials$gibbs_qualified[calib_trials$rep_id <= 100] &
          calib_trials$stan_qualified[calib_trials$rep_id <= 100]
      )
    )
  )
) %>%
  mutate(qualified_pct = 100 * n_qualified / n_trials)

write_export <- function(x, filename) {
  write.csv(x, file.path(export_dir, filename), row.names = FALSE, na = "")
}
write_export(study_design, "01_study_design.csv")
write_export(trial_censoring_summary, "02_censoring_summary.csv")
write_export(runtime_summary, "04_runtime_and_efficiency.csv")
write_export(runtime_long, "04_runtime_trial_level.csv")
write_export(speedup_summary, "04_speedup_summary.csv")
write_export(agreement_summary_table, "05_sampler_agreement.csv")
write_export(global_operating, "06_global_operating_characteristics.csv")
write_export(cohort_operating, "07_mixed_cohort_operating_characteristics.csv")
write_export(hetvar_groups, "14_hetvar_variance_groups.csv")
write_export(extreme_discrepancies, "15_extreme_discrepancies.csv")
write_export(paired, "20_paired_qualification.csv")
write_export(ess_index_summary, "25_ess_normalized_index.csv")
write_export(calib_qualification, "31_tcga_calibrated_qualification.csv")

cat("Wrote article evidence exports to ", export_dir, ".\n", sep = "")
