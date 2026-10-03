# Simulation study: data generation.
#
# Generates the deterministic trial datasets for every scenario in
# scenario_specs and writes data/sim_<scenario>_rep<nnn>.rds. Existing files are
# skipped, so the script resumes where it stopped.
#
# Usage: SIM_SCENARIOS=mixed SIM_N_REPS=500 Rscript 01_generate_data.R
# (on Windows PowerShell, set the two variables with $env: first).

dir.create("data", showWarnings = FALSE, recursive = TRUE)

#' Generate synthetic time-to-event data for a basket trial.
#'
#' @param n Total sample size across all baskets.
#' @param K Number of baskets.
#' @param scenario Clinical scenario identifier.
#' @param target_cens Nominal censoring proportion for the historical design.
#' @param seed Optional random seed for exact replicate reproduction.
#' @param error_dist Error distribution for event times ("normal", "gumbel", or "t").
#' @param t_df Degrees of freedom for Student-t errors.
#' @param censoring Censoring mechanism ("independent" or "informative").
#' @param cens_rho Correlation between censoring and event time on the error scale.
#' @param cens_sd Residual SD of the log censoring time.
#' @param cens_intercept Calibrated intercept for informative censoring; computed
#'   in a fixed-seed pilot when NULL.
#' @param sigma_by_basket Optional basket-specific residual SDs (length K).
#' @param group_sizes Optional exact basket sizes (length K; sum must equal n).
#' @param target_cens_by_basket Optional basket-specific censoring targets.
#'
#' @return A list containing the generated data, true parameters, and metadata.
generate_basket_data <- function(
  n = 135L,
  K = 9L,
  scenario = "mixed",
  target_cens = 0.30,
  seed = NULL,
  error_dist = c("normal", "gumbel", "t"),
  t_df = 5L,
  censoring = c("independent", "informative"),
  cens_rho = 0.5,
  cens_sd = 0.8,
  cens_intercept = NULL,
  sigma_by_basket = NULL,
  group_sizes = NULL,
  target_cens_by_basket = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (
    length(target_cens) != 1L ||
      !is.finite(target_cens) ||
      target_cens <= 0 ||
      target_cens >= 1
  ) {
    stop("target_cens must be a probability strictly between 0 and 1.")
  }
  if (
    !is.null(target_cens_by_basket) &&
      (length(target_cens_by_basket) != K ||
        any(!is.finite(target_cens_by_basket)) ||
        any(target_cens_by_basket <= 0) ||
        any(target_cens_by_basket >= 1))
  ) {
    stop("target_cens_by_basket must contain K probabilities in (0, 1).")
  }

  error_dist <- match.arg(error_dist)
  censoring <- match.arg(censoring)
  if (scenario == "weibull_misspec") {
    error_dist <- "gumbel"
  }
  if (scenario == "t_misspec") {
    error_dist <- "t"
  }
  if (scenario == "infcens_misspec") {
    censoring <- "informative"
  }
  if (scenario == "hetvar_misspec") {
    if (is.null(sigma_by_basket)) {
      sigma_by_basket <- seq(0.8, 1.6, length.out = K)
    }
    if (
      length(sigma_by_basket) != K ||
        any(!is.finite(sigma_by_basket)) ||
        any(sigma_by_basket <= 0)
    ) {
      stop("sigma_by_basket must contain K positive finite values.")
    }
  }

  sigma_true <- 1.2
  beta_true <- c(0.8, -0.5)
  theta_true <- numeric(K)
  tau_true <- 0.15

  if (scenario == "tcga_calib") {
    # Calibrated to the TCGA PanCancer Atlas application: fixed basket effects
    # equal to the exnexSurv posterior means, one active and one null covariate,
    # the application's residual SD, and per-basket sizes and censoring targets.
    n <- 4146L
    K <- 9L
    sigma_true <- 1.461
    beta_true <- c(-0.413, 0)
    tau_true <- 0.805
    theta_true <- c(
      4.287,
      3.984,
      5.425,
      4.329,
      4.804,
      3.559,
      3.962,
      4.581,
      4.233
    )
    mu_true <- mean(theta_true)
    if (is.null(group_sizes)) {
      group_sizes <- c(90L, 409L, 1071L, 284L, 568L, 182L, 522L, 510L, 510L)
    }
    if (is.null(target_cens_by_basket)) {
      target_cens_by_basket <- c(
        0.644,
        0.560,
        0.859,
        0.761,
        0.796,
        0.582,
        0.580,
        0.667,
        0.755
      )
    }
  } else if (scenario == "homogeneous") {
    mu_true <- 1.5
    tau_true <- 0.08
    theta_true <- rnorm(K, mean = mu_true, sd = tau_true)
  } else if (scenario == "mixed_divlarge") {
    # Mixed-efficacy effects with the three discordant baskets placed among the
    # largest arms, so that discordance is no longer confounded with basket
    # size. The random-number draws are identical to those of "mixed", which
    # lets the two scenarios be paired replicate by replicate.
    mu_true <- 1.5
    tau_true <- 0.15
    theta_true[1:3] <- c(-0.4, -0.6, -0.5)
    theta_true[4:9] <- rnorm(6L, mean = mu_true, sd = tau_true)
  } else if (
    scenario %in%
      c(
        "mixed",
        "weibull_misspec",
        "hetvar_misspec",
        "t_misspec",
        "infcens_misspec"
      )
  ) {
    mu_true <- 1.5
    tau_true <- 0.15
    theta_true[1:6] <- rnorm(6L, mean = mu_true, sd = tau_true)
    theta_true[7:9] <- c(-0.4, -0.6, -0.5)
  } else if (scenario == "heterogeneous") {
    mu_true <- 0.5
    tau_true <- 1.0
    theta_true <- seq(-1.2, 1.8, length.out = K)
  } else {
    stop("Unknown scenario: ", scenario)
  }

  if (
    !is.null(group_sizes) &&
      (length(group_sizes) != K ||
        any(!is.finite(group_sizes)) ||
        any(group_sizes != as.integer(group_sizes)) ||
        sum(group_sizes) != n ||
        any(group_sizes <= 0))
  ) {
    stop("group_sizes must contain K positive integers summing to n.")
  }

  if (!is.null(group_sizes)) {
    group <- sample(rep(seq_len(K), times = group_sizes))
  } else {
    basket_weights <- c(0.22, 0.18, 0.15, 0.12, 0.10, 0.08, 0.06, 0.05, 0.04)
    repeat {
      group <- sample(
        seq_len(K),
        size = n,
        replace = TRUE,
        prob = basket_weights
      )
      if (length(unique(group)) == K) {
        break
      }
    }
  }

  X <- matrix(rnorm(n * 2L), ncol = 2L)
  colnames(X) <- c("cov1", "cov2")
  mu_i <- theta_true[group] + drop(X %*% beta_true)
  if (error_dist == "gumbel") {
    sigma_use <- sigma_true * sqrt(6 / pi^2)
    log_T <- mu_i + sigma_use * (-log(-log(runif(n))))
  } else if (error_dist == "t") {
    sigma_use <- sigma_true * sqrt((t_df - 2) / t_df)
    log_T <- mu_i + sigma_use * rt(n, df = t_df)
  } else if (!is.null(sigma_by_basket)) {
    log_T <- mu_i + sigma_by_basket[group] * rnorm(n)
  } else {
    log_T <- rnorm(n, mean = mu_i, sd = sigma_true)
  }
  T_obs <- exp(log_T)

  if (censoring == "informative") {
    if (is.null(cens_intercept)) {
      if (exists("infcens_intercept_cache", envir = .GlobalEnv)) {
        cens_intercept <- get("infcens_intercept_cache", envir = .GlobalEnv)
      } else {
        cens_intercept <- calibrate_informative_censoring(
          # Use representative Mixed-efficacy locations so the calibrated
          # intercept is not anchored to one replicate's random draw.
          theta_true = c(rep(1.5, 6), -0.4, -0.6, -0.5),
          beta_true = beta_true,
          sigma_true = sigma_true,
          t_df = t_df,
          n = n,
          K = K,
          cens_rho = cens_rho,
          cens_sd = cens_sd
        )
        assign("infcens_intercept_cache", cens_intercept, envir = .GlobalEnv)
      }
    }
    log_C <- cens_intercept + cens_rho * (log_T - mu_i) + rnorm(n, 0, cens_sd)
    C_times <- exp(log_C)
  } else if (!is.null(target_cens_by_basket)) {
    rates <- calibrated_basket_rates(
      theta_true = theta_true,
      beta_true = beta_true,
      sigma_true = sigma_true,
      target_cens_by_basket = target_cens_by_basket,
      K = K
    )
    C_times <- rexp(n, rate = rates[group])
  } else {
    censoring_rate <- (1 - target_cens) / mean(T_obs)
    C_times <- rexp(n, rate = censoring_rate)
  }
  time <- pmin(T_obs, C_times)
  event <- as.numeric(T_obs <= C_times)

  data_df <- data.frame(
    time = time,
    event = event,
    group = group,
    cov1 = X[, 1L],
    cov2 = X[, 2L]
  )

  list(
    data = data_df,
    true_params = list(
      theta = theta_true,
      beta = beta_true,
      sigma2 = sigma_true^2,
      mu = mu_true,
      tau = tau_true
    ),
    metadata = list(
      scenario = scenario,
      n = n,
      K = K,
      target_cens = target_cens,
      censoring_mechanism = censoring,
      cens_rho = if (censoring == "informative") {
        cens_rho
      } else {
        NA_real_
      },
      cens_sd = if (censoring == "informative") {
        cens_sd
      } else {
        NA_real_
      },
      cens_intercept = if (censoring == "informative") {
        cens_intercept
      } else {
        NA_real_
      },
      target_cens_by_basket = if (is.null(target_cens_by_basket)) {
        rep(NA_real_, K)
      } else {
        target_cens_by_basket
      },
      group_sizes = as.numeric(table(factor(
        group,
        levels = seq_len(
          K
        )
      ))),
      observed_cens = mean(event == 0),
      error_distribution = error_dist,
      t_df = if (error_dist == "t") {
        t_df
      } else {
        NA_integer_
      },
      sigma_by_basket = if (is.null(sigma_by_basket)) {
        rep(sigma_true, K)
      } else {
        sigma_by_basket
      }
    )
  )
}

#' Calibrate the informative-censoring intercept to a target censoring rate.
#'
#' Uses a fixed-seed pilot so the calibrated intercept is reproducible.
calibrate_informative_censoring <- function(
  theta_true,
  beta_true,
  sigma_true,
  t_df = 5L,
  n = 135L,
  K = 9L,
  cens_rho = 0.5,
  cens_sd = 0.8,
  target = 0.30,
  tol = 0.02,
  pilot_reps = 60L,
  seed = 20260819L
) {
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  set.seed(seed)
  on.exit(
    {
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  basket_weights <- c(0.22, 0.18, 0.15, 0.12, 0.10, 0.08, 0.06, 0.05, 0.04)
  censoring_prop <- function(b) {
    props <- numeric(pilot_reps)
    for (r in seq_len(pilot_reps)) {
      group <- sample(
        seq_len(K),
        size = n,
        replace = TRUE,
        prob = basket_weights
      )
      if (length(unique(group)) < K) {
        group <- sample(
          seq_len(K),
          size = n,
          replace = TRUE,
          prob = basket_weights
        )
      }
      X <- matrix(rnorm(n * 2L), ncol = 2L)
      mu_i <- theta_true[group] + drop(X %*% beta_true)
      log_T <- rnorm(n, mean = mu_i, sd = sigma_true)
      log_C <- b + cens_rho * (log_T - mu_i) + rnorm(n, 0, cens_sd)
      props[r] <- mean(log_C < log_T)
    }
    mean(props)
  }
  lo <- -4
  hi <- 4
  for (i in seq_len(40L)) {
    mid <- (lo + hi) / 2
    p_mid <- censoring_prop(mid)
    if (abs(p_mid - target) <= tol) {
      break
    }
    if (p_mid > target) {
      lo <- mid
    } else {
      hi <- mid
    }
  }
  (lo + hi) / 2
}

#' Calibrate per-basket exponential censoring rates so the realized censoring
#' matches basket-specific targets under a heavy-tailed log-normal event-time
#' distribution. Uses a fixed-seed pilot and caches the result globally.
calibrated_basket_rates <- function(
  theta_true,
  beta_true,
  sigma_true,
  target_cens_by_basket,
  K = 9L,
  n_pilot = 5000L,
  seed = 20260819L
) {
  cache_name <- "tcga_calib_rates_cache"
  if (exists(cache_name, envir = .GlobalEnv, inherits = FALSE)) {
    return(get(cache_name, envir = .GlobalEnv, inherits = FALSE))
  }
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  set.seed(seed)
  on.exit(
    {
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  rates <- vapply(
    seq_len(K),
    function(k) {
      X <- matrix(rnorm(n_pilot * 2L), ncol = 2L)
      mu <- theta_true[k] + drop(X %*% beta_true)
      log_T <- rnorm(n_pilot, mean = mu, sd = sigma_true)
      T <- exp(log_T)
      target_c <- target_cens_by_basket[k]
      target_event <- 1 - target_c
      lo <- 1e-6
      hi <- 10
      for (i in seq_len(60L)) {
        mid <- sqrt(lo * hi)
        p <- mean(exp(-mid * T))
        if (abs(p - target_event) < 1e-3) {
          break
        }
        if (p > target_event) {
          lo <- mid
        } else {
          hi <- mid
        }
      }
      sqrt(lo * hi)
    },
    numeric(1)
  )
  assign(cache_name, rates, envir = .GlobalEnv)
  rates
}

scenario_specs <- c(
  homogeneous = 2000L,
  mixed = 2000L,
  heterogeneous = 2000L,
  weibull_misspec = 1000L,
  hetvar_misspec = 1000L,
  t_misspec = 1000L,
  infcens_misspec = 1000L,
  tcga_calib = 1000L,
  # Divergent-position robustness check (see 08_divergent_position.R); fitted
  # with the 'mixed' replicate seeds so the two scenarios can be paired.
  mixed_divlarge = 240L
)

# The extreme-value (max-Gumbel) misspecification scenario keeps the historical
# key weibull_misspec in data filenames and seeds; all public-facing labels use
# the corrected distributional name "Gumbel (extreme-value) misspecification".

scenarios <- names(scenario_specs)
n_reps_override <- Sys.getenv("SIM_N_REPS")
scenario_override <- Sys.getenv("SIM_SCENARIOS")
base_seed <- 1L

if (scenario_override != "") {
  requested_scenarios <- trimws(strsplit(scenario_override, ",", fixed = TRUE)[[
    1L
  ]])
  unknown_scenarios <- setdiff(requested_scenarios, scenarios)
  if (length(unknown_scenarios) > 0L) {
    stop("Unknown scenarios: ", paste(unknown_scenarios, collapse = ", "))
  }
  scenarios <- unique(requested_scenarios)
}

# Seed indices are pinned to their historical values so that removing the
# global-null design from the archive does not change any other scenario's
# random-number stream. mixed_divlarge reuses the 'mixed' index by design, so
# the two scenarios share their stream replicate by replicate.
seed_index_map <- c(
  homogeneous = 1L,
  mixed = 2L,
  heterogeneous = 3L,
  # global_null = 4L,
  weibull_misspec = 5L,
  hetvar_misspec = 6L,
  t_misspec = 7L,
  infcens_misspec = 8L,
  tcga_calib = 9L
)

for (scenario in scenarios) {
  n_reps <- if (n_reps_override != "") {
    as.integer(n_reps_override)
  } else {
    unname(scenario_specs[[scenario]])
  }
  for (replicate in seq_len(n_reps)) {
    seed_index <- if (scenario == "mixed_divlarge") {
      seed_index_map[["mixed"]]
    } else {
      seed_index_map[[scenario]]
    }
    replicate_seed <- base_seed + seed_index * 10000L + replicate
    output_file <- sprintf(
      "data/sim_%s_rep%03d.rds",
      scenario,
      replicate
    )
    if (file.exists(output_file)) {
      next
    }
    simulation <- generate_basket_data(
      scenario = scenario,
      seed = replicate_seed
    )
    saveRDS(simulation, file = output_file)
  }
}
