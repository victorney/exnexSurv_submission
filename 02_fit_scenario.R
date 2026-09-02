# Simulation study: fit all models for one scenario (made to run via Rscript).
# Install exnexSurv before running this script.

suppressPackageStartupMessages({
  library(survival)
  library(rstan)
  library(dplyr)
  library(tidyr)
  library(posterior)
  library(exnexSurv)
})

dir.create("fitted_models/scenarios", showWarnings = FALSE, recursive = TRUE)
rstan_options(auto_write = TRUE)

args <- commandArgs(trailingOnly = TRUE)
scenarios <- c(
  "homogeneous",
  "mixed",
  "heterogeneous",
  "global_null",
  "weibull_misspec",
  "hetvar_misspec",
  "t_misspec",
  "infcens_misspec"
)

if (length(args) != 1L || !args[[1L]] %in% scenarios) {
  stop(
    "Provide one scenario: ",
    paste(scenarios, collapse = ", "),
    ". Example: Rscript 02_fit_scenario.R mixed"
  )
}

scen <- args[[1L]]
scenario_index <- match(scen, scenarios)

scenario_n_reps <- c(
  homogeneous = 2000L,
  mixed = 2000L,
  heterogeneous = 2000L,
  global_null = 2000L,
  weibull_misspec = 1000L,
  hetvar_misspec = 1000L,
  t_misspec = 1000L,
  infcens_misspec = 1000L
)

n_reps_override <- Sys.getenv("FIT_N_REPS")

n_reps <- if (n_reps_override != "") {
  as.integer(n_reps_override)
} else {
  unname(scenario_n_reps[[scen]])
}

n_iter <- 2000
n_warmup <- 1000
n_chains <- 4
fit_seed_base <- 1000000L + scenario_index * 10000L

options(mc.cores = 4)

# Compile the comparison models once per scenario worker.
stan_exnex <- stan_model("models/exnex_marginal.stan")
stan_pooled <- stan_model("models/pooled_marginal.stan")
stan_unpooled <- stan_model("models/unpooled_marginal.stan")

# Posterior summary helpers
extract_stan_summary <- function(fit, param_name, K_rep) {
  na_df <- data.frame(
    mean = rep(NA_real_, K_rep),
    sd = rep(NA_real_, K_rep),
    q025 = rep(NA_real_, K_rep),
    q975 = rep(NA_real_, K_rep),
    rhat = rep(NA_real_, K_rep),
    ess = rep(NA_real_, K_rep)
  )

  if (is.null(fit)) {
    return(na_df)
  }

  summ <- tryCatch(
    summary(fit, pars = param_name)$summary,
    error = function(e) NULL
  )

  if (is.null(summ) || nrow(summ) == 0) {
    return(na_df)
  }

  col_orig <- colnames(summ)
  if (nrow(summ) == 1 && K_rep > 1) {
    summ <- matrix(rep(summ, K_rep), nrow = K_rep, byrow = TRUE)
    colnames(summ) <- col_orig
  }

  ess_col <- if ("n_eff" %in% col_orig) {
    "n_eff"
  } else if ("Bulk_ESS" %in% col_orig) {
    "Bulk_ESS"
  } else {
    col_orig[1]
  }

  data.frame(
    mean = as.numeric(summ[, "mean"]),
    sd = as.numeric(summ[, "sd"]),
    q025 = as.numeric(summ[, "2.5%"]),
    q975 = as.numeric(summ[, "97.5%"]),
    rhat = as.numeric(summ[, "Rhat"]),
    ess = as.numeric(summ[, ess_col]),
    row.names = NULL
  )
}

extract_exnex_summary <- function(fit, K_rep, n_chains, n_post_per_chain) {
  draws_mat <- fit$draws
  theta_cols <- paste0("theta_", 1:K_rep)
  res_list <- vector("list", K_rep)

  for (k in 1:K_rep) {
    vec_draws <- draws_mat[, theta_cols[k]]
    mat_chains <- matrix(vec_draws, nrow = n_post_per_chain, ncol = n_chains)
    res_list[[k]] <- data.frame(
      mean = mean(vec_draws),
      sd = sd(vec_draws),
      q025 = as.numeric(quantile(vec_draws, probs = 0.025)),
      q975 = as.numeric(quantile(vec_draws, probs = 0.975)),
      rhat = tryCatch(
        as.numeric(posterior::rhat(mat_chains)),
        error = function(e) NA_real_
      ),
      ess = tryCatch(
        as.numeric(posterior::ess_bulk(mat_chains)),
        error = function(e) NA_real_
      )
    )
  }

  res_df <- do.call(rbind, res_list)
  res_df
}

# Resume from an existing scenario checkpoint.
scen_file <- sprintf("fitted_models/scenarios/results_scenario_%s.rds", scen)

if (file.exists(scen_file)) {
  scen_results_df <- readRDS(scen_file) %>%
    select(-any_of(c("grupo", "exnex_pexch")))

  scen_results_df$target_cens <- 0.30

  scen_results_df$observed_cens <- 1 -
    scen_results_df$n_events / scen_results_df$n_total

  if (any(!is.finite(scen_results_df$observed_cens))) {
    stop("Checkpoint contains non-finite censoring proportions.")
  }

  completed_reps <- unique(scen_results_df$rep_id)
  scen_results <- split(scen_results_df, scen_results_df$rep_id)
} else {
  scen_results <- list()
  completed_reps <- c()
}

for (r in seq_len(n_reps)) {
  if (r %in% completed_reps) {
    next
  }

  data_file <- sprintf("data/sim_%s_rep%03d.rds", scen, r)
  if (!file.exists(data_file)) {
    stop(sprintf("Dataset for '%s' replicate %d not found.", scen, r))
  }

  sim <- readRDS(data_file)
  df <- sim$data
  observed_cens <- mean(df$event == 0)
  K <- sim$metadata$K
  rep_seed <- fit_seed_base + r * 10L

  basket_counts <- data.frame(group = 1:K) %>%
    left_join(
      df %>%
        group_by(group) %>%
        summarise(n_total = n(), n_events = sum(event == 1), .groups = "drop"),
      by = "group"
    ) %>%
    mutate(
      n_total = tidyr::replace_na(n_total, 0L),
      n_events = tidyr::replace_na(n_events, 0L)
    ) %>%
    arrange(group)

  # 1. exnexSurv (Gibbs DA)
  t_exnex_val <- NA_real_
  fit_exnex <- tryCatch(
    {
      t_start <- proc.time()[[3L]]
      f <- exnex_surv(
        Surv(time, event) ~ group + cov1 + cov2,
        data = df,
        iter = n_iter,
        warmup = n_warmup,
        chains = n_chains,
        parallel_chains = n_chains,
        seed = rep_seed + 1L
      )
      t_exnex_val <- proc.time()[[3L]] - t_start
      f
    },
    error = function(e) NULL
  )

  # 2. Stan Models
  df_obs <- df[df$event == 1, ]
  df_cens <- df[df$event == 0, ]
  stan_data <- list(
    N_obs = nrow(df_obs),
    N_cens = nrow(df_cens),
    K = K,
    P = 2,
    log_y_obs = log(df_obs$time),
    log_y_cens = array(log(df_cens$time), dim = nrow(df_cens)),
    group_obs = df_obs$group,
    group_cens = array(df_cens$group, dim = nrow(df_cens)),
    X_obs = as.matrix(df_obs[, c("cov1", "cov2")]),
    X_cens = as.matrix(df_cens[, c("cov1", "cov2")])
  )

  t_stan_exnex_val <- NA_real_
  fit_stan_exnex <- tryCatch(
    {
      t_start <- proc.time()[[3L]]
      f <- sampling(
        stan_exnex,
        data = stan_data,
        iter = n_iter,
        warmup = n_warmup,
        chains = n_chains,
        seed = rep_seed + 2L,
        refresh = 0,
        show_messages = FALSE,
        open_progress = FALSE
      )
      t_stan_exnex_val <- proc.time()[[3L]] - t_start
      f
    },
    error = function(e) NULL
  )

  t_stan_pooled_val <- NA_real_
  fit_stan_pooled <- tryCatch(
    {
      t_start <- proc.time()[[3L]]
      f <- sampling(
        stan_pooled,
        data = stan_data,
        iter = n_iter,
        warmup = n_warmup,
        chains = n_chains,
        seed = rep_seed + 3L,
        refresh = 0,
        show_messages = FALSE,
        open_progress = FALSE
      )
      t_stan_pooled_val <- proc.time()[[3L]] - t_start
      f
    },
    error = function(e) NULL
  )

  t_stan_unp_val <- NA_real_
  fit_stan_unpooled <- tryCatch(
    {
      t_start <- proc.time()[[3L]]
      f <- sampling(
        stan_unpooled,
        data = stan_data,
        iter = n_iter,
        warmup = n_warmup,
        chains = n_chains,
        seed = rep_seed + 4L,
        refresh = 0,
        show_messages = FALSE,
        open_progress = FALSE
      )
      t_stan_unp_val <- proc.time()[[3L]] - t_start
      f
    },
    error = function(e) NULL
  )

  n_post_chain <- n_iter - n_warmup
  exnex_summ <- if (!is.null(fit_exnex)) {
    extract_exnex_summary(fit_exnex, K, n_chains, n_post_chain)
  } else {
    data.frame(
      mean = rep(NA_real_, K),
      sd = rep(NA_real_, K),
      q025 = rep(NA_real_, K),
      q975 = rep(NA_real_, K),
      rhat = rep(NA_real_, K),
      ess = rep(NA_real_, K)
    )
  }
  stan_exnex_summ <- extract_stan_summary(fit_stan_exnex, "theta", K)
  stan_pool_summ <- extract_stan_summary(fit_stan_pooled, "theta", K)
  stan_unp_summ <- extract_stan_summary(fit_stan_unpooled, "theta", K)

  res_iter <- data.frame(
    scenario = scen,
    rep_id = r,
    target_cens = 0.30,
    observed_cens = observed_cens,
    group = 1:K,
    n_total = basket_counts$n_total,
    n_events = basket_counts$n_events,
    theta_true = sim$true_params$theta,
    time_exnex = t_exnex_val,
    time_stan_exnex = t_stan_exnex_val,
    time_stan_pooled = t_stan_pooled_val,
    time_stan_unp = t_stan_unp_val,
    exnex_mean = exnex_summ$mean,
    exnex_sd = exnex_summ$sd,
    exnex_q025 = exnex_summ$q025,
    exnex_q975 = exnex_summ$q975,
    exnex_rhat = exnex_summ$rhat,
    exnex_ess = exnex_summ$ess,
    stan_exnex.mean = stan_exnex_summ$mean,
    stan_exnex.sd = stan_exnex_summ$sd,
    stan_exnex.q025 = stan_exnex_summ$q025,
    stan_exnex.q975 = stan_exnex_summ$q975,
    stan_exnex.rhat = stan_exnex_summ$rhat,
    stan_exnex.ess = stan_exnex_summ$ess,
    stan_pooled.mean = stan_pool_summ$mean,
    stan_pooled.sd = stan_pool_summ$sd,
    stan_pooled.q025 = stan_pool_summ$q025,
    stan_pooled.q975 = stan_pool_summ$q975,
    stan_pooled.rhat = stan_pool_summ$rhat,
    stan_pooled.ess = stan_pool_summ$ess,
    stan_unpooled.mean = stan_unp_summ$mean,
    stan_unpooled.sd = stan_unp_summ$sd,
    stan_unpooled.q025 = stan_unp_summ$q025,
    stan_unpooled.q975 = stan_unp_summ$q975,
    stan_unpooled.rhat = stan_unp_summ$rhat,
    stan_unpooled.ess = stan_unp_summ$ess,
    row.names = NULL
  )

  scen_results[[as.character(r)]] <- res_iter

  if (length(scen_results) %% 25 == 0 || r == n_reps) {
    checkpoint <- do.call(rbind, scen_results)
    saveRDS(checkpoint, scen_file)
    rm(sim, df, fit_exnex, fit_stan_exnex, fit_stan_pooled, fit_stan_unpooled)
    gc(verbose = FALSE)
  }
}
