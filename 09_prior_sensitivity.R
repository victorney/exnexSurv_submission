# Prior and borrowing-control sensitivity for exnexSurv.
#
# Fits the Mixed Efficacy design under nine prior configurations on the same
# 500 archived replicates and writes the article export 12_prior_sensitivity.csv.
# Deterministic: the data are read from data/sim_mixed_rep*.rds and each fit
# uses seed = 70000 + 10000 * configuration index + replicate.

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(posterior)
  library(exnexSurv)
})

dir.create("fitted_models/scenarios", showWarnings = FALSE, recursive = TRUE)

n_reps <- 500
n_iter <- 2000
n_warmup <- 1000
n_chains <- 4

configs <- list(
  baseline = list(),
  p_mix_03 = list(p_mix = 0.3),
  p_mix_07 = list(p_mix = 0.7),
  v_nex_10 = list(v_nex = 10),
  v_nex_100 = list(v_nex = 100),
  tau_ig33 = list(a_tau = 3, b_tau = 3),
  sigma_ig33 = list(a_sigma = 3, b_sigma = 3),
  mu_prior = list(m_mu = 0, v_mu = 100),
  borrowing_strong = list(p_mix = 0.7, v_nex = 10)
)
config_labels <- c(
  baseline = "Default priors",
  p_mix_03 = "p_exch = 0.3",
  p_mix_07 = "p_exch = 0.7",
  v_nex_10 = "v_nex = 10",
  v_nex_100 = "v_nex = 100",
  tau_ig33 = "tau2 ~ IG(3,3)",
  sigma_ig33 = "sigma2 ~ IG(3,3)",
  mu_prior = "mu ~ N(0, 100)",
  borrowing_strong = "p_exch = 0.7, v_nex = 10"
)

extract_summary <- function(fit, K, n_chains, n_post_per_chain) {
  if (is.null(fit)) {
    return(data.frame(
      mean = rep(NA_real_, K),
      sd = rep(NA_real_, K),
      q025 = rep(NA_real_, K),
      q975 = rep(NA_real_, K),
      rhat = rep(NA_real_, K),
      ess = rep(NA_real_, K)
    ))
  }
  draws_mat <- fit$draws
  theta_cols <- paste0("theta_", seq_len(K))
  do.call(
    rbind,
    lapply(seq_len(K), function(k) {
      vec <- draws_mat[, theta_cols[k]]
      mat <- matrix(vec, nrow = n_post_per_chain, ncol = n_chains)
      data.frame(
        mean = mean(vec),
        sd = sd(vec),
        q025 = as.numeric(quantile(vec, probs = 0.025)),
        q975 = as.numeric(quantile(vec, probs = 0.975)),
        rhat = tryCatch(as.numeric(posterior::rhat(mat)), error = function(e) {
          NA_real_
        }),
        ess = tryCatch(
          as.numeric(posterior::ess_bulk(mat)),
          error = function(e) NA_real_
        )
      )
    })
  )
}

for (cfg_name in names(configs)) {
  out_file <- sprintf(
    "fitted_models/scenarios/results_prior_sensitivity_%s.rds",
    cfg_name
  )
  if (file.exists(out_file)) {
    scen_results_df <- readRDS(out_file)
    completed_reps <- unique(scen_results_df$rep_id)
    scen_results <- split(scen_results_df, scen_results_df$rep_id)
  } else {
    scen_results <- list()
    completed_reps <- c()
  }

  cfg_idx <- match(cfg_name, names(configs))
  for (r in seq_len(n_reps)) {
    if (r %in% completed_reps) {
      next
    }
    data_file <- sprintf("data/sim_mixed_rep%03d.rds", r)
    if (!file.exists(data_file)) {
      stop(
        "Missing ",
        data_file,
        ". Generate the Mixed-efficacy replicates first: ",
        "SIM_SCENARIOS=mixed SIM_N_REPS=",
        n_reps,
        " Rscript 01_generate_data.R"
      )
    }
    sim <- readRDS(data_file)
    df <- sim$data
    K <- sim$metadata$K

    fit <- tryCatch(
      exnex_surv(
        Surv(time, event) ~ group + cov1 + cov2,
        data = df,
        priors = configs[[cfg_name]],
        iter = n_iter,
        warmup = n_warmup,
        chains = n_chains,
        # exnexSurv 1.3.x takes a chain count; 1.4.0 takes a logical flag.
        parallel_chains = if (utils::packageVersion("exnexSurv") >= "1.4.0") {
          TRUE
        } else {
          n_chains
        },
        seed = 70000L + cfg_idx * 10000L + r
      ),
      error = function(e) NULL
    )
    summ <- extract_summary(fit, K, n_chains, n_iter - n_warmup)

    res_iter <- data.frame(
      config = cfg_name,
      rep_id = r,
      group = seq_len(K),
      theta_true = sim$true_params$theta,
      n_total = as.numeric(table(factor(df$group, levels = seq_len(K)))),
      n_events = as.numeric(table(factor(
        df$group[df$event == 1],
        levels = seq_len(K)
      ))),
      mean = summ$mean,
      sd = summ$sd,
      q025 = summ$q025,
      q975 = summ$q975,
      rhat = summ$rhat,
      ess = summ$ess
    )
    scen_results[[as.character(r)]] <- res_iter

    if (length(scen_results) %% 25 == 0 || r == n_reps) {
      checkpoint <- do.call(rbind, scen_results)
      saveRDS(checkpoint, out_file)
      rm(sim, df, fit)
      gc(verbose = FALSE)
    }
  }
}

results <- bind_rows(lapply(
  names(configs),
  function(cfg_name) {
    readRDS(sprintf(
      "fitted_models/scenarios/results_prior_sensitivity_%s.rds",
      cfg_name
    ))
  }
))

results <- results %>%
  mutate(
    config = factor(
      config,
      levels = names(configs),
      labels = unname(config_labels)
    ),
    error = mean - theta_true,
    interval_width = q975 - q025,
    covered_95 = theta_true >= q025 & theta_true <= q975,
    finite_summary = is.finite(mean) & is.finite(q025) & is.finite(q975),
    ordered_interval = is.finite(q025) & is.finite(q975) & q025 <= q975,
    rhat_lt_1_01 = is.finite(rhat) & rhat < 1.01,
    ess_ge_400 = is.finite(ess) & ess >= 400,
    qualified = finite_summary & ordered_interval & rhat_lt_1_01 & ess_ge_400,
    cohort = case_when(
      group <= 6L ~ "Responsive 1-6",
      group >= 7L ~ "Resistant 7-9",
      TRUE ~ NA_character_
    )
  )

trial_status <- results %>%
  group_by(config, rep_id) %>%
  summarise(
    n_baskets = n_distinct(group),
    qualified_trial = n_baskets == 9L & all(qualified),
    .groups = "drop"
  )

results <- results %>%
  left_join(trial_status, by = c("config", "rep_id"))

cohort_summary <- results %>%
  group_by(config, cohort) %>%
  summarise(
    n_basket_evaluations = n(),
    bias = mean(error, na.rm = TRUE),
    error_sd = sd(error, na.rm = TRUE),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    bias2_plus_error_sd2 = bias^2 + error_sd^2,
    coverage_95 = mean(covered_95, na.rm = TRUE) * 100,
    mean_interval_width = mean(interval_width, na.rm = TRUE),
    median_abs_error = median(abs(error), na.rm = TRUE),
    p90_abs_error = as.numeric(quantile(abs(error), 0.90, na.rm = TRUE)),
    p95_abs_error = as.numeric(quantile(abs(error), 0.95, na.rm = TRUE)),
    max_abs_error = max(abs(error), na.rm = TRUE),
    prop_abs_error_gt_0_5 = mean(abs(error) > 0.5, na.rm = TRUE) * 100,
    prop_abs_error_gt_1 = mean(abs(error) > 1, na.rm = TRUE) * 100,
    qualified_basket_rate = mean(qualified, na.rm = TRUE) * 100,
    .groups = "drop"
  )

qualification_summary <- trial_status %>%
  group_by(config) %>%
  summarise(
    n_trials = n(),
    qualified_trial_rate = mean(qualified_trial) * 100,
    .groups = "drop"
  )

decision_metrics <- results %>%
  group_by(config, cohort) %>%
  summarise(
    n_confident = sum((q025 > 0) | (q975 < 0), na.rm = TRUE),
    sign_error_rate = mean(
      (mean > 0) != (theta_true > 0),
      na.rm = TRUE
    ) *
      100,
    confident_sign_error = ifelse(
      n_confident > 0,
      sum(
        ((q025 > 0) & (theta_true <= 0)) |
          ((q975 < 0) & (theta_true > 0)),
        na.rm = TRUE
      ) /
        n_confident *
        100,
      NA_real_
    ),
    .groups = "drop"
  )

export <- cohort_summary %>%
  left_join(qualification_summary, by = "config") %>%
  left_join(decision_metrics, by = c("config", "cohort")) %>%
  arrange(config, cohort)

dir.create(
  "fitted_models/article_exports",
  showWarnings = FALSE,
  recursive = TRUE
)
write.csv(
  export,
  "fitted_models/article_exports/12_prior_sensitivity.csv",
  row.names = FALSE,
  na = ""
)

print(
  as.data.frame(export) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))),
  row.names = FALSE
)

invisible(NULL)
