# Benchmark: centered versus non-centered Stan EXNEX parameterizations.
# Fits both parameterizations on a fixed subset of archived replicates and
# compares runtime, bulk ESS per second, divergences, and maximum treedepth.

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
subset_size <- as.integer(Sys.getenv("NC_SUBSET_SIZE", unset = "100"))
options(mc.cores = 4)

stan_centered <- stan_model("models/exnex_marginal.stan")
stan_noncentered <- stan_model("models/exnex_marginal_noncentered.stan")

extract_theta_ess <- function(fit, K) {
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

extract_theta_means <- function(fit, K) {
  summ <- summary(fit, pars = "theta")$summary
  as.numeric(summ[, "mean"])
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
    error = function(e) NULL
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

scenarios <- c(mixed = "mixed", heterogeneous = "heterogeneous")
base_seed <- 90000L
results_list <- list()

for (scenario_index in seq_along(scenarios)) {
  scen <- unname(scenarios[[scenario_index]])
  for (r in seq_len(subset_size)) {
    data_file <- sprintf("data/sim_%s_rep%03d.rds", scen, r)
    if (!file.exists(data_file)) {
      warning("Missing dataset: ", data_file)
      next
    }
    sim <- readRDS(data_file)
    stan_data <- build_stan_data(sim)
    rep_seed <- base_seed + scenario_index * 10000L + r

    centered <- fit_trial(stan_centered, stan_data, rep_seed)
    noncentered <- fit_trial(stan_noncentered, stan_data, rep_seed + 5000L)

    K <- sim$metadata$K
    row_out <- data.frame(
      scenario = scen,
      rep_id = r,
      variant = c("centered", "noncentered"),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        runtime_seconds = c(centered$seconds, noncentered$seconds),
        mean_theta_ess = c(
          if (!is.null(centered$fit)) {
            extract_theta_ess(centered$fit, K)
          } else {
            NA_real_
          },
          if (!is.null(noncentered$fit)) {
            extract_theta_ess(noncentered$fit, K)
          } else {
            NA_real_
          }
        ),
        ess_per_second = mean_theta_ess / runtime_seconds,
        n_divergent = c(
          if (!is.null(centered$fit)) {
            sampler_diagnostics(centered$fit)$n_divergent
          } else {
            NA_real_
          },
          if (!is.null(noncentered$fit)) {
            sampler_diagnostics(noncentered$fit)$n_divergent
          } else {
            NA_real_
          }
        ),
        max_treedepth = c(
          if (!is.null(centered$fit)) {
            sampler_diagnostics(centered$fit)$max_treedepth
          } else {
            NA_real_
          },
          if (!is.null(noncentered$fit)) {
            sampler_diagnostics(noncentered$fit)$max_treedepth
          } else {
            NA_real_
          }
        )
      )
    results_list[[length(results_list) + 1L]] <- row_out

    if (r %% 10 == 0 || r == subset_size) {
      cat(sprintf(
        "noncentered benchmark: %s rep %d/%d done\n",
        scen,
        r,
        subset_size
      ))
    }
  }
}

results <- bind_rows(results_list)

agreement <- results %>%
  select(scenario, rep_id, variant, mean_theta_ess, runtime_seconds) %>%
  pivot_wider(
    id_cols = c(scenario, rep_id),
    names_from = variant,
    values_from = c(mean_theta_ess, runtime_seconds)
  )

write.csv(
  results,
  "fitted_models/article_exports/18_noncentered_stan.csv",
  row.names = FALSE,
  na = ""
)
write.csv(
  agreement,
  "fitted_models/article_exports/18_noncentered_stan_paired.csv",
  row.names = FALSE,
  na = ""
)

summary_table <- results %>%
  group_by(scenario, variant) %>%
  summarise(
    n_trials = n(),
    runtime_mean = mean(runtime_seconds, na.rm = TRUE),
    ess_per_second_mean = mean(ess_per_second, na.rm = TRUE),
    ess_per_second_median = median(ess_per_second, na.rm = TRUE),
    total_divergences = sum(n_divergent, na.rm = TRUE),
    max_treedepth = max(max_treedepth, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Non-centered Stan benchmark ===\n")
print(as.data.frame(summary_table), row.names = FALSE)
