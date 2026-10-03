# Benchmark: Stan EXNEX geometry on the Complete Heterogeneity stress subset.
# Fits the centered parameterization on the first 100 archived heterogeneous
# replicates, records runtime, effective sample size, divergent transitions
# and maximum treedepth, and derives the Gibbs-versus-Stan stress summary
# stratified by divergence burden (export 22).
#
# Inputs: data/sim_heterogeneous_rep001..100.rds (generate them first with
#   SIM_SCENARIOS=heterogeneous SIM_N_REPS=100 Rscript 01_generate_data.R)
#   and fitted_models/scenarios/results_scenario_heterogeneous.rds.
# Outputs: fitted_models/article_exports/18_stan_stress_benchmark.csv and
#   fitted_models/article_exports/22_gibbs_on_stan_stress.csv.

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(tidyr)
})

dir.create(
  "fitted_models/article_exports",
  showWarnings = FALSE,
  recursive = TRUE
)
rstan_options(auto_write = TRUE)

n_iter <- 2000
n_warmup <- 1000
n_chains <- 4
subset_size <- as.integer(Sys.getenv("STRESS_SUBSET_SIZE", unset = "100"))
# The heterogeneous scenario occupied the second position of the original
# scenario vector, so its fit seeds keep that historical offset.
base_seed <- 90000L + 20000L
options(mc.cores = 4)

stan_centered <- stan_model("models/exnex_marginal.stan")

extract_theta_ess <- function(fit) {
  summ <- summary(fit, pars = "theta")$summary
  ess_col <- if ("n_eff" %in% colnames(summ)) {
    "n_eff"
  } else if ("Bulk_ESS" %in% colnames(summ)) {
    "Bulk_ESS"
  } else {
    colnames(summ)[1]
  }
  mean(as.numeric(summ[, ess_col]))
}

sampler_diagnostics <- function(fit) {
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  divergent <- vapply(
    sampler_params,
    function(x) sum(x[, "divergent__"]),
    numeric(1)
  )
  treedepth <- vapply(
    sampler_params,
    function(x) max(x[, "treedepth__"]),
    numeric(1)
  )
  list(n_divergent = sum(divergent), max_treedepth = max(treedepth))
}

fit_trial <- function(model, stan_data, seed) {
  t0 <- proc.time()[[3L]]
  fit <- tryCatch(
    sampling(
      model,
      data = stan_data,
      iter = n_iter,
      warmup = n_warmup,
      chains = n_chains,
      seed = seed,
      refresh = 0,
      show_messages = FALSE,
      open_progress = FALSE
    ),
    error = function(e) {
      warning("Stan fit failed for seed ", seed, ": ", conditionMessage(e))
      NULL
    }
  )
  list(fit = fit, seconds = proc.time()[[3L]] - t0)
}

build_stan_data <- function(sim) {
  df <- sim$data
  df_obs <- df[df$event == 1, ]
  df_cens <- df[df$event == 0, ]
  list(
    N_obs = nrow(df_obs),
    N_cens = nrow(df_cens),
    K = sim$metadata$K,
    P = 2,
    log_y_obs = log(df_obs$time),
    log_y_cens = array(log(df_cens$time), dim = nrow(df_cens)),
    group_obs = df_obs$group,
    group_cens = array(df_cens$group, dim = nrow(df_cens)),
    X_obs = as.matrix(df_obs[, c("cov1", "cov2")]),
    X_cens = as.matrix(df_cens[, c("cov1", "cov2")])
  )
}

results_list <- list()
for (r in seq_len(subset_size)) {
  data_file <- sprintf("data/sim_heterogeneous_rep%03d.rds", r)
  if (!file.exists(data_file)) {
    stop(
      "Missing ",
      data_file,
      ". Generate the stress subset first: ",
      "SIM_SCENARIOS=heterogeneous SIM_N_REPS=",
      subset_size,
      " Rscript 01_generate_data.R"
    )
  }
  sim <- readRDS(data_file)
  stan_data <- build_stan_data(sim)
  fit <- fit_trial(stan_centered, stan_data, base_seed + r)

  results_list[[length(results_list) + 1L]] <- data.frame(
    scenario = "heterogeneous",
    rep_id = r,
    runtime_seconds = fit$seconds,
    mean_theta_ess = if (!is.null(fit$fit)) {
      extract_theta_ess(fit$fit)
    } else {
      NA_real_
    },
    n_divergent = if (!is.null(fit$fit)) {
      sampler_diagnostics(fit$fit)$n_divergent
    } else {
      NA_real_
    },
    max_treedepth = if (!is.null(fit$fit)) {
      sampler_diagnostics(fit$fit)$max_treedepth
    } else {
      NA_real_
    }
  )

  if (r %% 10 == 0 || r == subset_size) {
    cat(sprintf("stress benchmark: rep %d/%d done\n", r, subset_size))
  }
}

benchmark <- bind_rows(results_list) %>%
  mutate(ess_per_second = mean_theta_ess / runtime_seconds)

write.csv(
  benchmark,
  "fitted_models/article_exports/18_stan_stress_benchmark.csv",
  row.names = FALSE,
  na = ""
)

# Stress summary stratified by the number of divergent transitions in the
# centered fit, using the archived heterogeneous fits for the Gibbs and Stan
# qualification and error summaries.
ckpt <- readRDS("fitted_models/scenarios/results_scenario_heterogeneous.rds")

basket_long <- bind_rows(
  ckpt %>%
    transmute(
      rep_id,
      group,
      theta_true,
      method_id = "gibbs_exnex",
      estimate = exnex_mean,
      q025 = exnex_q025,
      q975 = exnex_q975,
      rhat = exnex_rhat,
      ess = exnex_ess
    ),
  ckpt %>%
    transmute(
      rep_id,
      group,
      theta_true,
      method_id = "stan_exnex",
      estimate = stan_exnex.mean,
      q025 = stan_exnex.q025,
      q975 = stan_exnex.q975,
      rhat = stan_exnex.rhat,
      ess = stan_exnex.ess
    )
) %>%
  mutate(
    finite = is.finite(estimate) & is.finite(q025) & is.finite(q975),
    ordered = is.finite(q025) & is.finite(q975) & q025 <= q975,
    qualified = finite &
      ordered &
      is.finite(rhat) &
      rhat < 1.01 &
      is.finite(ess) &
      ess >= 400,
    error = estimate - theta_true
  )

wide_paired <- basket_long %>%
  select(rep_id, group, method_id, estimate, error, rhat, ess, qualified) %>%
  pivot_wider(
    names_from = method_id,
    values_from = c(estimate, error, rhat, ess, qualified),
    names_glue = "{.value}_{method_id}"
  ) %>%
  filter(!is.na(estimate_gibbs_exnex), !is.na(estimate_stan_exnex))

trial_metrics <- wide_paired %>%
  group_by(rep_id) %>%
  summarise(
    gibbs_qualified = all(qualified_gibbs_exnex),
    gibbs_max_rhat = max(rhat_gibbs_exnex),
    gibbs_min_ess = min(ess_gibbs_exnex),
    gibbs_median_abs_error = median(abs(error_gibbs_exnex)),
    gibbs_trial_rmse = sqrt(mean(error_gibbs_exnex^2)),
    stan_qualified = all(qualified_stan_exnex),
    stan_trial_rmse = sqrt(mean(error_stan_exnex^2)),
    .groups = "drop"
  )

stress_summary <- trial_metrics %>%
  inner_join(benchmark %>% select(rep_id, n_divergent), by = "rep_id") %>%
  mutate(
    divergence_stratum = cut(
      n_divergent,
      breaks = c(-Inf, 0, 10, 100, Inf),
      labels = c("0", "1-10", "11-100", ">100"),
      right = TRUE
    )
  ) %>%
  group_by(divergence_stratum) %>%
  summarise(
    scenario = "heterogeneous",
    n_trials = n(),
    stan_divergences_total = sum(n_divergent),
    gibbs_qualified_pct = 100 * mean(gibbs_qualified, na.rm = TRUE),
    stan_qualified_pct = 100 * mean(stan_qualified, na.rm = TRUE),
    gibbs_median_trial_rmse = median(gibbs_trial_rmse, na.rm = TRUE),
    stan_median_trial_rmse = median(stan_trial_rmse, na.rm = TRUE),
    gibbs_median_max_rhat = median(gibbs_max_rhat, na.rm = TRUE),
    gibbs_median_min_ess = median(gibbs_min_ess, na.rm = TRUE),
    gibbs_median_abs_error = median(gibbs_median_abs_error, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  select(scenario, divergence_stratum, everything())

write.csv(
  stress_summary,
  "fitted_models/article_exports/22_gibbs_on_stan_stress.csv",
  row.names = FALSE,
  na = ""
)

cat("\n=== Stan stress benchmark ===\n")
print(
  benchmark %>%
    summarise(
      n_trials = n(),
      runtime_mean = mean(runtime_seconds, na.rm = TRUE),
      total_divergences = sum(n_divergent, na.rm = TRUE),
      max_treedepth = max(max_treedepth, na.rm = TRUE)
    ) %>%
    as.data.frame(),
  row.names = FALSE
)
cat("\n=== Gibbs on Stan stress subset ===\n")
print(as.data.frame(stress_summary), row.names = FALSE)
