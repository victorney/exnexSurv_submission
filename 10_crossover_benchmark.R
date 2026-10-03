# Crossover benchmark separating basket size from censoring.
# Eight baskets: the largest grows from 400 to 1,200 patients while the other
# seven stay at 371, 343, 314, 286, 257, 229 and 200. Five exchangeable
# effects are drawn around 1.5 and three discordant effects (-0.4, -0.6, -0.5)
# sit in baskets 6-8. Two censoring regimes are compared: 30% in every basket
# and a TCGA-like regime (86% in the largest basket, 65% in the others, 71%
# overall). The script fits exnexSurv and the marginalized centered Stan model
# on three datasets per size and writes the accuracy-matched index ratios
# behind Table 4 of the article.
#
# Inputs: none beyond models/exnex_marginal.stan.
# Outputs: fitted_models/article_exports/29_crossover_benchmark.csv (per
#   dataset) and 30_crossover_benchmark_summary.csv (ratios by size/regime).
# Environment: CROSSOVER_ITER (default 2000 total iterations),
#   CROSSOVER_DATASETS (default 3).

suppressPackageStartupMessages({
  library(survival)
  library(rstan)
  library(exnexSurv)
  library(posterior)
  library(dplyr)
  library(tidyr)
})

dir.create(
  "fitted_models/article_exports",
  showWarnings = FALSE,
  recursive = TRUE
)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

n_iter <- as.integer(Sys.getenv("CROSSOVER_ITER", unset = "2000"))
n_warmup <- n_iter %/% 2L
n_chains <- 4L
n_datasets <- as.integer(Sys.getenv("CROSSOVER_DATASETS", unset = "3"))

sizes_other <- c(371L, 343L, 314L, 286L, 257L, 229L, 200L)
size_grid <- c(400L, 600L, 800L, 1000L, 1200L)
discordant <- c(-0.4, -0.6, -0.5)
beta_true <- c(0.8, -0.5)
sigma_true <- 1.2
censoring_regimes <- list(
  moderate = rep(0.30, 8L),
  tcga_like = c(0.86, rep(0.65, 7L))
)

gen_data <- function(seed, largest, cens_targets) {
  set.seed(seed)
  sizes <- c(largest, sizes_other)
  n <- sum(sizes)
  group <- rep(seq_along(sizes), times = sizes)
  cov1 <- rnorm(n)
  cov2 <- rnorm(n)
  theta <- c(rnorm(5L, 1.5, 0.15), discordant)
  eta <- theta[group] + beta_true[1] * cov1 + beta_true[2] * cov2
  event_time <- exp(rnorm(n, eta, sigma_true))
  cens_time <- numeric(n)
  for (b in seq_along(sizes)) {
    idx <- group == b
    lambda <- uniroot(
      function(l) mean(exp(-l * event_time[idx])) - (1 - cens_targets[b]),
      c(1e-8, 5)
    )$root
    cens_time[idx] <- rexp(sum(idx), lambda)
  }
  data.frame(
    time = pmin(event_time, cens_time),
    event = as.integer(event_time <= cens_time),
    group = group,
    cov1 = cov1,
    cov2 = cov2
  )
}

stan_data_from <- function(df) {
  obs <- df[df$event == 1, ]
  cens <- df[df$event == 0, ]
  list(
    N_obs = nrow(obs),
    N_cens = nrow(cens),
    K = 8L,
    P = 2L,
    log_y_obs = log(obs$time),
    log_y_cens = array(log(cens$time), dim = nrow(cens)),
    group_obs = obs$group,
    group_cens = array(cens$group, dim = nrow(cens)),
    X_obs = as.matrix(obs[, c("cov1", "cov2")]),
    X_cens = as.matrix(cens[, c("cov1", "cov2")])
  )
}

fit_exnex <- function(df, seed) {
  exnex_surv(
    Surv(time, event) ~ group + cov1 + cov2,
    data = df,
    iter = n_iter,
    warmup = n_warmup,
    chains = n_chains,
    # exnexSurv 1.3.x takes the number of chains to run concurrently;
    # version 1.4.0 replaced that with a logical flag.
    parallel_chains = if (utils::packageVersion("exnexSurv") >= "1.4.0") {
      TRUE
    } else {
      n_chains
    },
    seed = seed
  )
}

fit_stan <- function(stan_data, seed) {
  sampling(
    stan_mod,
    data = stan_data,
    iter = n_iter,
    warmup = n_warmup,
    chains = n_chains,
    cores = n_chains,
    seed = seed,
    refresh = 0,
    show_messages = FALSE,
    open_progress = FALSE
  )
}

gibbs_metrics <- function(fit) {
  draws <- fit$draws
  n_post <- nrow(draws) / n_chains
  ess <- vapply(
    seq_len(8L),
    function(k) {
      vec <- draws[[paste0("theta_", k)]]
      posterior::ess_bulk(matrix(vec, nrow = n_post, ncol = n_chains))
    },
    numeric(1)
  )
  list(min_ess = min(ess), mean_ess = mean(ess))
}

stan_metrics <- function(fit) {
  summ <- summary(fit, pars = "theta")$summary
  ess <- as.numeric(summ[, "n_eff"])
  sp <- get_sampler_params(fit, inc_warmup = FALSE)
  list(
    min_ess = min(ess),
    mean_ess = mean(ess),
    n_divergent = sum(vapply(
      sp,
      function(x) sum(x[, "divergent__"]),
      numeric(1)
    ))
  )
}

stan_mod <- stan_model("models/exnex_marginal.stan")

results <- list()
for (regime in names(censoring_regimes)) {
  regime_index <- match(regime, names(censoring_regimes))
  for (size_index in seq_along(size_grid)) {
    for (dataset in seq_len(n_datasets)) {
      data_seed <- if (regime == "moderate") {
        7000L + 100L * size_index + dataset
      } else {
        5000L + 100L * size_index + dataset
      }
      df <- gen_data(
        data_seed,
        size_grid[size_index],
        censoring_regimes[[regime]]
      )
      stan_data <- stan_data_from(df)

      t0 <- proc.time()[[3L]]
      gibbs_fit <- tryCatch(
        fit_exnex(
          df,
          seed = 100000L + 10000L * regime_index + 100L * size_index + dataset
        ),
        error = function(e) {
          warning("exnexSurv failed: ", conditionMessage(e))
          NULL
        }
      )
      gibbs_seconds <- proc.time()[[3L]] - t0

      t0 <- proc.time()[[3L]]
      stan_fit <- tryCatch(
        fit_stan(
          stan_data,
          seed = 200000L + 10000L * regime_index + 100L * size_index + dataset
        ),
        error = function(e) {
          warning("Stan fit failed: ", conditionMessage(e))
          NULL
        }
      )
      stan_seconds <- proc.time()[[3L]] - t0

      gibbs <- if (!is.null(gibbs_fit)) {
        gibbs_metrics(gibbs_fit)
      } else {
        list(min_ess = NA_real_, mean_ess = NA_real_)
      }
      stan <- if (!is.null(stan_fit)) {
        stan_metrics(stan_fit)
      } else {
        list(min_ess = NA_real_, mean_ess = NA_real_, n_divergent = NA_real_)
      }

      results[[length(results) + 1L]] <- data.frame(
        regime = regime,
        largest_basket = size_grid[size_index],
        dataset = dataset,
        n_total = nrow(df),
        observed_censoring = 1 - mean(df$event),
        method = c("gibbs", "stan"),
        runtime_seconds = c(gibbs_seconds, stan_seconds),
        min_ess = c(gibbs$min_ess, stan$min_ess),
        mean_ess = c(gibbs$mean_ess, stan$mean_ess),
        n_divergent = c(NA_real_, stan$n_divergent)
      )
      cat(sprintf(
        "%s L=%d ds=%d | gibbs %.1fs | stan %.1fs\n",
        regime,
        size_grid[size_index],
        dataset,
        gibbs_seconds,
        stan_seconds
      ))
    }
  }
}

benchmark <- bind_rows(results) %>%
  mutate(index_t_star = runtime_seconds * 400 / min_ess)

summary_table <- benchmark %>%
  group_by(regime, largest_basket, method) %>%
  summarise(
    index_t_star = mean(index_t_star, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = method,
    values_from = index_t_star,
    names_prefix = "index_"
  ) %>%
  mutate(ratio_stan_over_gibbs = index_stan / index_gibbs) %>%
  arrange(regime, largest_basket)

write.csv(
  benchmark,
  "fitted_models/article_exports/29_crossover_benchmark.csv",
  row.names = FALSE,
  na = ""
)
write.csv(
  summary_table,
  "fitted_models/article_exports/30_crossover_benchmark_summary.csv",
  row.names = FALSE,
  na = ""
)

cat("\n=== Crossover index ratios (Stan / Gibbs) ===\n")
print(as.data.frame(summary_table), row.names = FALSE)
