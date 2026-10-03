# Divergent-position check for the Mixed Efficacy design.
#
# Mixed Efficacy puts the three discordant effects in the three smallest baskets
# (expected sizes 8, 7, 5); the mixed_divlarge scenario puts the same three
# effects in the three largest baskets (expected 30, 24, 20). Both scenarios are
# generated from the same replicate seeds, so paired trials share their
# covariates, their group assignments and their random draws, and the only
# difference between them is where the discordance sits.
#
# Reproduce with:
#   SIM_SCENARIOS=mixed_divlarge Rscript 01_generate_data.R
#   FIT_REP_FROM=1   FIT_REP_TO=120 FIT_SHARD=1 Rscript 02_fit_scenario.R mixed_divlarge
#   FIT_REP_FROM=121 FIT_REP_TO=240 FIT_SHARD=2 Rscript 02_fit_scenario.R mixed_divlarge
#   Rscript 08_divergent_position.R      # writes exports 26, 27 and 28
#
# Populations follow 03_analyze_article.R: jointly qualified trials for the two
# EXNEX implementations, with the paired set defined as the replicates that
# qualified in both designs.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(purrr)
})

export_dir <- "fitted_models/article_exports"
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

read_scenario <- function(scen) {
  files <- list.files(
    "fitted_models/scenarios",
    pattern = sprintf("^results_scenario_%s(_shard[0-9]+)?\\.rds$", scen),
    full.names = TRUE
  )
  if (length(files) == 0L) {
    stop("No checkpoint found for scenario '", scen, "'.")
  }
  bind_rows(lapply(files, readRDS)) %>%
    distinct(scenario, rep_id, group, .keep_all = TRUE)
}

methods <- tribble(
  ~method_id        , ~method                , ~estimate         , ~q025             , ~q975 , ~rhat , ~ess ,
  "gibbs_exnex"     , "exnexSurv (Gibbs DA)" , "exnex_mean"      , "exnex_q025"      ,
  "exnex_q975"      , "exnex_rhat"           , "exnex_ess"       ,
  "stan_exnex"      , "Stan EXNEX (NUTS)"    , "stan_exnex.mean" , "stan_exnex.q025" ,
  "stan_exnex.q975" , "stan_exnex.rhat"      , "stan_exnex.ess"
)

cohort_labels <- function(df) {
  case_when(
    df$scenario == "mixed" & df$group <= 6L ~ "Exchangeable baskets",
    df$scenario == "mixed" & df$group >= 7L ~ "Discordant baskets",
    df$scenario == "mixed_divlarge" & df$group <= 3L ~ "Discordant baskets",
    df$scenario == "mixed_divlarge" & df$group >= 4L ~ "Exchangeable baskets"
  )
}

as_long <- function(df) {
  bind_rows(lapply(seq_len(nrow(methods)), function(i) {
    m <- methods[i, ]
    tibble(
      scenario = df$scenario,
      rep_id = df$rep_id,
      group = df$group,
      cohort = cohort_labels(df),
      n_total = df$n_total,
      n_events = df$n_events,
      theta_true = df$theta_true,
      estimate = df[[m$estimate]],
      q025 = df[[m$q025]],
      q975 = df[[m$q975]],
      rhat = df[[m$rhat]],
      ess = df[[m$ess]],
      method_id = m$method_id,
      method = m$method
    )
  }))
}

scenario_labels <- c(
  mixed = "Mixed Efficacy",
  mixed_divlarge = "Divergent-position"
)
scenario_order <- c("mixed", "mixed_divlarge")

long <- bind_rows(
  as_long(read_scenario("mixed")),
  as_long(read_scenario("mixed_divlarge"))
) %>%
  mutate(
    scenario_label = unname(scenario_labels[scenario]),
    finite = is.finite(estimate) & is.finite(q025) & is.finite(q975),
    ordered = is.finite(q025) & is.finite(q975) & q025 <= q975,
    qualified_basket = finite &
      ordered &
      is.finite(rhat) &
      rhat < 1.01 &
      is.finite(ess) &
      ess >= 400
  )

# Trial-level qualification, mirroring 03_analyze_article.R.
trial_status <- long %>%
  group_by(scenario, rep_id, method_id) %>%
  summarise(
    n_baskets = n_distinct(group),
    qualified_trial = n_baskets == 9L & all(qualified_basket),
    .groups = "drop"
  )

joint_trials <- trial_status %>%
  group_by(scenario, rep_id) %>%
  summarise(
    jointly_qualified = n() == nrow(methods) & all(qualified_trial),
    .groups = "drop"
  )

long <- long %>%
  left_join(trial_status, by = c("scenario", "rep_id", "method_id")) %>%
  left_join(joint_trials, by = c("scenario", "rep_id")) %>%
  filter(finite) %>%
  mutate(
    population = case_when(
      jointly_qualified ~ "jointly_qualified",
      qualified_trial ~ "qualified_trials",
      TRUE ~ "all_fits"
    ),
    error = estimate - theta_true,
    covered = theta_true >= q025 & theta_true <= q975,
    width = q975 - q025
  )

# Common evaluation set: replicates in which both implementations qualified in
# both designs, so every reported row rests on exactly the same trials.
paired_reps <- joint_trials %>%
  pivot_wider(names_from = scenario, values_from = jointly_qualified) %>%
  filter(!is.na(mixed), !is.na(mixed_divlarge), mixed, mixed_divlarge) %>%
  pull(rep_id)

paired_long <- long %>%
  filter(rep_id %in% paired_reps, population == "jointly_qualified")

trial_metrics <- paired_long %>%
  group_by(scenario, scenario_label, method_id, method, cohort, rep_id) %>%
  summarise(
    trial_bias = mean(error),
    trial_mse = mean(error^2),
    trial_coverage = mean(covered),
    trial_width = mean(width),
    .groups = "drop"
  )

cohort_summary <- paired_long %>%
  group_by(scenario, scenario_label, method_id, method, cohort) %>%
  summarise(
    mean_basket_size = mean(n_total),
    mean_events = mean(n_events),
    n_basket_trials = n(),
    bias = mean(error),
    coverage_95 = 100 * mean(covered),
    mean_interval_width = mean(width),
    .groups = "drop"
  ) %>%
  left_join(
    trial_metrics %>%
      group_by(scenario, method_id, cohort) %>%
      summarise(
        n_trials = n_distinct(rep_id),
        mcse_bias = sd(trial_bias) / sqrt(n()),
        rmse = sqrt(mean(trial_mse)),
        mcse_rmse = sd(trial_mse) / sqrt(n()) / (2 * sqrt(mean(trial_mse))),
        mcse_coverage = 100 * sd(trial_coverage) / sqrt(n()),
        mcse_interval_width = sd(trial_width) / sqrt(n()),
        .groups = "drop"
      ),
    by = c("scenario", "method_id", "cohort")
  ) %>%
  mutate(across(c(bias, rmse, coverage_95, mean_interval_width), \(x) {
    round(x, 4)
  }))

# Paired contrast inside the discordant cohort: same replicates, same draws.
paired <- trial_metrics %>%
  filter(cohort == "Discordant baskets") %>%
  select(
    method_id,
    method,
    rep_id,
    scenario,
    trial_bias,
    trial_mse,
    trial_coverage
  ) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(trial_bias, trial_mse, trial_coverage)
  ) %>%
  filter(!is.na(trial_bias_mixed), !is.na(trial_bias_mixed_divlarge)) %>%
  mutate(
    delta_bias = trial_bias_mixed_divlarge - trial_bias_mixed,
    delta_coverage = 100 *
      (trial_coverage_mixed_divlarge - trial_coverage_mixed),
    delta_rmse = sqrt(trial_mse_mixed_divlarge) - sqrt(trial_mse_mixed)
  )

paired_summary <- paired %>%
  group_by(method_id, method) %>%
  summarise(
    n_pairs = n(),
    mean_delta_bias = mean(delta_bias),
    mcse_delta_bias = sd(delta_bias) / sqrt(n()),
    mean_delta_rmse = mean(delta_rmse),
    mcse_delta_rmse = sd(delta_rmse) / sqrt(n()),
    mean_delta_coverage = mean(delta_coverage),
    mcse_delta_coverage = sd(delta_coverage) / sqrt(n()),
    p_delta_bias = t.test(delta_bias)$p.value,
    p_delta_coverage = t.test(delta_coverage)$p.value,
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric) & !starts_with("p_"), \(x) round(x, 5)))

# Qualification rates over every replicate available in both designs, which is
# the number the article quotes (the paired set above is conditioned on it).
common_reps <- long %>%
  group_by(scenario) %>%
  summarise(reps = list(unique(rep_id)), .groups = "drop") %>%
  pull(reps) %>%
  purrr::reduce(intersect)

qualification <- trial_status %>%
  filter(rep_id %in% common_reps) %>%
  group_by(scenario, method_id) %>%
  summarise(
    n_trials = n(),
    qualified_pct = 100 * mean(qualified_trial),
    .groups = "drop"
  ) %>%
  mutate(
    scenario_label = unname(scenario_labels[scenario]),
    method = unname(setNames(methods$method, methods$method_id)[method_id])
  )

write_csv(cohort_summary, file.path(export_dir, "26_divergent_position.csv"))
write_csv(
  paired_summary,
  file.path(export_dir, "27_divergent_position_paired.csv")
)
write_csv(
  qualification,
  file.path(export_dir, "28_divergent_position_qualification.csv")
)

cat("\n== Paired replicates:", length(paired_reps), "==\n")
cat("\n== Cohort operating characteristics (jointly qualified) ==\n")
print(
  cohort_summary %>%
    transmute(
      design = scenario_label,
      cohort,
      method,
      size = round(mean_basket_size, 1),
      events = round(mean_events, 1),
      bias,
      rmse,
      coverage = coverage_95,
      mcse_cov = round(mcse_coverage, 2),
      n = n_trials
    ) %>%
    as.data.frame()
)

cat(
  "\n== Paired discordant-cohort contrast (divergent-position minus Mixed Efficacy) ==\n"
)
print(
  paired_summary %>%
    transmute(
      method,
      n_pairs,
      d_bias = mean_delta_bias,
      se_bias = mcse_delta_bias,
      d_rmse = mean_delta_rmse,
      d_cov = mean_delta_coverage,
      se_cov = mcse_delta_coverage,
      p_bias = signif(p_delta_bias, 2),
      p_cov = signif(p_delta_coverage, 2)
    ) %>%
    as.data.frame()
)

cat("\n== Trial qualification rates on the same replicates ==\n")
print(as.data.frame(qualification))
