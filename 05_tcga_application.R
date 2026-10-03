# TCGA PanCancer Atlas application: EXNEX log-normal AFT models for basket trials.
# Downloads open-access clinical data from the cBioPortal datahub, builds
# overall-survival baskets by cancer type, and fits the four comparison models.

suppressPackageStartupMessages({
  library(survival)
  library(rstan)
  library(dplyr)
  library(tidyr)
  library(posterior)
  library(exnexSurv)
})

data_dir <- "data/tcga"
raw_dir <- file.path(data_dir, "raw")
export_dir <- "fitted_models/tcga_application"
article_export_dir <- "fitted_models/article_exports"
figure_dir <- "figures"

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(article_export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

studies <- c(
  "acc",
  "blca",
  "brca",
  "cesc",
  "coadread",
  "esca",
  "hnsc",
  "kirc",
  "lgg"
)
datahub_commit <- "0cc9138746c08b304f8dac92c31983e0ef44af1d"
base_url <- paste0(
  "https://media.githubusercontent.com/media/cBioPortal/datahub/",
  datahub_commit,
  "/public"
)
expected_md5 <- c(
  acc = "d381cb8b44c1493c5dfe535869c08ee5",
  blca = "a0750c2867c29b81ccd2fc0eda575bd1",
  brca = "373427b6cbd2fb16215e9c4d37e6ca6f",
  cesc = "21ff9491741b931bde8e225456fdc906",
  coadread = "eb695f1d1361d738ae51c8d329f51a8c",
  esca = "3128f014231c42243668875daaa1b279",
  hnsc = "716005db3070a8435580d1740132879c",
  kirc = "99ced25c11cfc1c1e5e68995c7ada28c",
  lgg = "67ae327aeb596753f8017f76193ecb3c"
)

n_iter <- 2000L
n_warmup <- 1000L
n_chains <- 4L
seed <- 42L

rstan_options(auto_write = TRUE)
options(mc.cores = n_chains)

basket_levels <- c(
  "ACC",
  "BLCA",
  "BRCA",
  "CESC",
  "COADREAD",
  "ESCA",
  "HNSC",
  "KIRC",
  "LGG"
)

method_levels <- c(
  "Complete pooling",
  "No pooling (stratified)",
  "Stan EXNEX (NUTS)",
  "exnexSurv (Gibbs DA)"
)
download_clinical <- function(study) {
  dest <- file.path(
    raw_dir,
    paste0(study, "_tcga_pan_can_atlas_2018_clinical_patient.txt")
  )
  if (file.exists(dest) && file.info(dest)$size > 0) {
    if (!identical(unname(tools::md5sum(dest)), expected_md5[[study]])) {
      stop("Existing TCGA file has an unexpected checksum: ", dest)
    }
    return(dest)
  }
  url <- paste0(
    base_url,
    "/",
    study,
    "_tcga_pan_can_atlas_2018/data_clinical_patient.txt"
  )
  download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
  if (grepl("git-lfs", readLines(dest, n = 1L))) {
    stop("Downloaded a Git-LFS pointer; the media URL should return raw text.")
  }
  if (!identical(unname(tools::md5sum(dest)), expected_md5[[study]])) {
    stop("Downloaded TCGA file has an unexpected checksum: ", dest)
  }
  dest
}

read_clinical <- function(file) {
  read.delim(
    file,
    header = TRUE,
    sep = "\t",
    comment.char = "#",
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  ) %>%
    select(PATIENT_ID, CANCER_TYPE_ACRONYM, AGE, OS_STATUS, OS_MONTHS)
}

clinical_files <- vapply(studies, download_clinical, character(1))
clinical <- bind_rows(lapply(clinical_files, read_clinical))

tcga <- clinical %>%
  mutate(
    basket = ifelse(
      CANCER_TYPE_ACRONYM %in% c("COAD", "READ"),
      "COADREAD",
      CANCER_TYPE_ACRONYM
    ),
    event = as.integer(startsWith(OS_STATUS, "1")),
    os_months = as.numeric(OS_MONTHS),
    age = as.numeric(AGE)
  ) %>%
  filter(
    is.finite(os_months),
    os_months > 0,
    !is.na(event),
    event %in% c(0L, 1L),
    is.finite(age)
  ) %>%
  select(patient_id = PATIENT_ID, basket, os_months, event, age) %>%
  mutate(
    basket = factor(basket, levels = basket_levels),
    age_std = as.numeric(scale(age))
  )

K <- nlevels(tcga$basket)
P <- 1L

saveRDS(tcga, file.path(data_dir, "tcga_clinical.rds"))

basket_summary <- tcga %>%
  group_by(basket) %>%
  summarise(
    n = n(),
    n_events = sum(event),
    censoring = 1 - mean(event),
    median_age = median(age),
    median_followup = median(os_months),
    .groups = "drop"
  )

stan_data <- list(
  N_obs = sum(tcga$event == 1L),
  N_cens = sum(tcga$event == 0L),
  K = K,
  P = P,
  log_y_obs = log(tcga$os_months[tcga$event == 1L]),
  log_y_cens = log(tcga$os_months[tcga$event == 0L]),
  group_obs = as.integer(tcga$basket)[tcga$event == 1L],
  group_cens = as.integer(tcga$basket)[tcga$event == 0L],
  X_obs = matrix(tcga$age_std[tcga$event == 1L], ncol = P),
  X_cens = matrix(tcga$age_std[tcga$event == 0L], ncol = P)
)

fit_gibbs_exnex <- function() {
  t0 <- proc.time()[[3L]]
  fit <- tryCatch(
    exnex_surv(
      Surv(os_months, event) ~ basket + age_std,
      data = tcga,
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
    ),
    error = function(e) {
      warning("exnexSurv failed: ", conditionMessage(e))
      NULL
    }
  )
  list(fit = fit, seconds = proc.time()[[3L]] - t0)
}

fit_stan_model <- function(model) {
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
      warning("Stan fit failed: ", conditionMessage(e))
      NULL
    }
  )
  list(fit = fit, seconds = proc.time()[[3L]] - t0)
}

cat("Fitting exnexSurv (Gibbs DA)...\n")
gibbs <- fit_gibbs_exnex()

stan_exnex_model <- stan_model("models/exnex_marginal.stan")
stan_pooled_model <- stan_model("models/pooled_marginal.stan")
stan_unpooled_model <- stan_model("models/unpooled_marginal.stan")

cat("Fitting Stan EXNEX (NUTS)...\n")
stan_exnex <- fit_stan_model(stan_exnex_model)
cat("Fitting complete pooling...\n")
stan_pooled <- fit_stan_model(stan_pooled_model)
cat("Fitting no pooling...\n")
stan_unpooled <- fit_stan_model(stan_unpooled_model)

extract_stan_summary <- function(fit, param_name) {
  if (is.null(fit)) {
    return(NULL)
  }
  summ <- tryCatch(
    summary(fit, pars = param_name)$summary,
    error = function(e) NULL
  )
  if (is.null(summ) || nrow(summ) == 0) {
    return(NULL)
  }
  col_orig <- colnames(summ)
  ess_col <- if ("n_eff" %in% col_orig) {
    "n_eff"
  } else if ("Bulk_ESS" %in% col_orig) {
    "Bulk_ESS"
  } else {
    col_orig[1]
  }
  data.frame(
    estimate = as.numeric(summ[, "mean"]),
    posterior_sd = as.numeric(summ[, "sd"]),
    q025 = as.numeric(summ[, "2.5%"]),
    q975 = as.numeric(summ[, "97.5%"]),
    rhat = as.numeric(summ[, "Rhat"]),
    bulk_ess = as.numeric(summ[, ess_col]),
    row.names = NULL
  )
}

summarize_draws_df <- function(draws) {
  n_post <- nrow(draws) / n_chains
  if (n_post != as.integer(n_post)) {
    stop("Draws do not divide evenly into chains.")
  }
  bind_rows(lapply(colnames(draws), function(parameter) {
    values <- draws[[parameter]]
    chain_mat <- matrix(values, nrow = n_post, ncol = n_chains)
    data.frame(
      parameter = parameter,
      estimate = mean(values),
      posterior_sd = sd(values),
      q025 = as.numeric(quantile(values, 0.025)),
      q975 = as.numeric(quantile(values, 0.975)),
      rhat = tryCatch(
        as.numeric(posterior::rhat(chain_mat)),
        error = function(e) NA_real_
      ),
      bulk_ess = tryCatch(
        as.numeric(posterior::ess_bulk(chain_mat)),
        error = function(e) NA_real_
      )
    )
  }))
}

assemble_block <- function(method_id, method, parameter, basket_id, summ) {
  if (is.null(summ) || nrow(summ) == 0) {
    return(NULL)
  }
  data.frame(
    method_id = method_id,
    method = method,
    parameter = parameter,
    basket_id = basket_id,
    summ,
    stringsAsFactors = FALSE
  )
}

theta_block <- function(method_id, method, summ, K_rep) {
  if (nrow(summ) == 1L && K_rep > 1L) {
    summ <- summ[rep(1L, K_rep), , drop = FALSE]
  }
  if (nrow(summ) != K_rep) {
    stop("theta summary has ", nrow(summ), " rows; expected ", K_rep, ".")
  }
  assemble_block(method_id, method, "theta", seq_len(K_rep), summ)
}

gibbs_long <- NULL
if (!is.null(gibbs$fit)) {
  gibbs_draws <- summarize_draws_df(gibbs$fit$draws)
  gibbs_long <- bind_rows(
    theta_block(
      "gibbs_exnex",
      "exnexSurv (Gibbs DA)",
      gibbs_draws %>%
        filter(grepl("^theta_", parameter)) %>%
        select(-parameter),
      K
    ),
    assemble_block(
      "gibbs_exnex",
      "exnexSurv (Gibbs DA)",
      "beta",
      NA_integer_,
      gibbs_draws %>% filter(parameter == "beta_1") %>% select(-parameter)
    ),
    assemble_block(
      "gibbs_exnex",
      "exnexSurv (Gibbs DA)",
      "sigma2",
      NA_integer_,
      gibbs_draws %>% filter(parameter == "sigma2") %>% select(-parameter)
    )
  )
}

stan_long <- function(method_id, method, fit, with_exnex_hyper = FALSE) {
  if (is.null(fit)) {
    return(NULL)
  }
  bind_rows(
    theta_block(method_id, method, extract_stan_summary(fit, "theta"), K),
    assemble_block(
      method_id,
      method,
      "beta",
      NA_integer_,
      extract_stan_summary(fit, "beta")
    ),
    assemble_block(
      method_id,
      method,
      "sigma2",
      NA_integer_,
      extract_stan_summary(fit, "sigma2")
    ),
    if (with_exnex_hyper) {
      bind_rows(
        assemble_block(
          method_id,
          method,
          "mu",
          NA_integer_,
          extract_stan_summary(fit, "mu")
        ),
        assemble_block(
          method_id,
          method,
          "tau2",
          NA_integer_,
          extract_stan_summary(fit, "tau2")
        )
      )
    }
  )
}

posterior_summary <- bind_rows(
  gibbs_long,
  stan_long(
    "stan_exnex",
    "Stan EXNEX (NUTS)",
    stan_exnex$fit,
    with_exnex_hyper = TRUE
  ),
  stan_long("stan_pooled", "Complete pooling", stan_pooled$fit),
  stan_long("stan_unpooled", "No pooling (stratified)", stan_unpooled$fit)
)

if (is.null(posterior_summary) || nrow(posterior_summary) == 0) {
  stop("All model fits failed; no posterior summaries available.")
}

posterior_summary <- posterior_summary %>%
  mutate(
    method = factor(method, levels = method_levels),
    basket = ifelse(
      is.na(basket_id),
      NA_character_,
      basket_levels[basket_id]
    )
  )

diagnostics <- posterior_summary %>%
  select(method_id, method, parameter, basket, rhat, bulk_ess)

runtime <- data.frame(
  method_id = c("gibbs_exnex", "stan_exnex", "stan_pooled", "stan_unpooled"),
  method = method_levels[c(4, 3, 1, 2)],
  seconds = c(
    gibbs$seconds,
    stan_exnex$seconds,
    stan_pooled$seconds,
    stan_unpooled$seconds
  )
)

source_metadata <- data.frame(
  study = studies,
  url = paste0(
    base_url,
    "/",
    studies,
    "_tcga_pan_can_atlas_2018/data_clinical_patient.txt"
  ),
  file = clinical_files,
  rows = vapply(clinical_files, function(f) nrow(read_clinical(f)), integer(1)),
  stringsAsFactors = FALSE
)

source_metadata <- source_metadata %>%
  mutate(
    commit_sha = datahub_commit,
    file_checksum_md5 = unname(tools::md5sum(file))
  )

write_export <- function(x, filename) {
  write.csv(x, file.path(export_dir, filename), row.names = FALSE, na = "")
}

write_export(basket_summary, "01_basket_summary.csv")
write_export(posterior_summary, "02_posterior_summary.csv")
write_export(diagnostics, "03_mcmc_diagnostics.csv")
write_export(runtime, "04_runtime.csv")
write_export(source_metadata, "00_source_metadata.csv")

saveRDS(
  list(
    data = tcga,
    basket_summary = basket_summary,
    posterior = posterior_summary,
    diagnostics = diagnostics,
    runtime = runtime,
    gibbs = gibbs$fit,
    stan_exnex = stan_exnex$fit,
    stan_pooled = stan_pooled$fit,
    stan_unpooled = stan_unpooled$fit
  ),
  file.path(export_dir, "tcga_results.rds")
)

# Restricted mean survival time at 60 months from posterior draws (95% CrI).
normalize_theta_cols <- function(d) {
  nm <- names(d)
  theta_idx <- grepl("^theta", nm)
  nm[theta_idx] <- paste0("theta_", seq_len(sum(theta_idx)))
  names(d) <- nm
  d
}

rmst_summary <- function(draws_df, K, t_star = 60) {
  d <- normalize_theta_cols(as.data.frame(draws_df))
  theta_cols <- paste0("theta_", seq_len(K))
  sigma2_vec <- d$sigma2
  if (is.null(sigma2_vec) || !all(is.finite(sigma2_vec))) {
    stop("sigma2 draws are required for RMST computation.")
  }
  bind_rows(lapply(seq_len(K), function(k) {
    mu <- d[[theta_cols[k]]]
    sigma <- sqrt(sigma2_vec)
    a <- (log(t_star) - mu) / sigma
    rmst <- exp(mu + sigma^2 / 2) * pnorm(a - sigma) + t_star * pnorm(-a)
    data.frame(
      basket = basket_levels[k],
      rmst_mean = mean(rmst),
      rmst_q025 = as.numeric(quantile(rmst, 0.025)),
      rmst_q975 = as.numeric(quantile(rmst, 0.975)),
      stringsAsFactors = FALSE
    )
  }))
}

validate_rmst_formula <- function(draws_df, K, t_star = 60) {
  d <- normalize_theta_cols(as.data.frame(draws_df))
  check_draws <- seq(1, nrow(d), length.out = 20)
  for (k in seq_len(min(K, 3L))) {
    for (i in check_draws) {
      mu <- d[[paste0("theta_", k)]][i]
      sigma <- sqrt(d$sigma2[i])
      a <- (log(t_star) - mu) / sigma
      analytic <- exp(mu + sigma^2 / 2) * pnorm(a - sigma) + t_star * pnorm(-a)
      numeric_val <- integrate(
        function(u) 1 - pnorm((log(u) - mu) / sigma),
        lower = 0,
        upper = t_star
      )$value
      if (abs(analytic - numeric_val) > 0.5) {
        warning(
          "RMST analytic formula mismatch at basket ",
          k,
          ", draw ",
          i,
          ": ",
          analytic,
          " vs ",
          numeric_val
        )
      }
    }
  }
  invisible(TRUE)
}

rmst_rows <- list()
if (!is.null(gibbs$fit)) {
  validate_rmst_formula(gibbs$fit$draws, K)
  rmst_rows[["exnexSurv (Gibbs DA)"]] <- rmst_summary(gibbs$fit$draws, K) %>%
    mutate(method = "exnexSurv (Gibbs DA)")
}
if (!is.null(stan_exnex$fit)) {
  rmst_rows[["Stan EXNEX (NUTS)"]] <- rmst_summary(
    posterior::as_draws_df(stan_exnex$fit),
    K
  ) %>%
    mutate(method = "Stan EXNEX (NUTS)")
}
if (length(rmst_rows) > 0L) {
  tcga_rmst <- bind_rows(rmst_rows) %>%
    select(method, basket, rmst_mean, rmst_q025, rmst_q975)
  write_export(tcga_rmst, "16_tcga_rmst.csv")
  write.csv(
    tcga_rmst,
    file.path(article_export_dir, "16_tcga_rmst.csv"),
    row.names = FALSE,
    na = ""
  )
}

theta_summary <- posterior_summary %>%
  filter(parameter == "theta") %>%
  left_join(
    basket_summary %>% select(basket, n, n_events),
    by = "basket"
  )

basket_order <- theta_summary %>%
  filter(method_id == "gibbs_exnex") %>%
  arrange(estimate) %>%
  pull(basket)

theta_summary <- theta_summary %>%
  mutate(basket = factor(basket, levels = basket_order))

panel_a_data <- theta_summary %>% filter(method_id != "stan_pooled")
pooled_band <- theta_summary %>%
  filter(method_id == "stan_pooled") %>%
  slice(1)

shrinkage_data <- theta_summary %>%
  filter(method_id %in% c("gibbs_exnex", "stan_unpooled")) %>%
  select(basket, method_id, estimate, q025, q975) %>%
  pivot_wider(
    names_from = method_id,
    values_from = c(estimate, q025, q975)
  ) %>%
  left_join(basket_summary %>% select(basket, n_events), by = "basket")

beta_summary <- posterior_summary %>%
  filter(parameter == "beta")

sigma_summary <- posterior_summary %>%
  filter(parameter == "sigma2") %>%
  mutate(
    estimate = sqrt(estimate),
    q025 = sqrt(q025),
    q975 = sqrt(q975)
  )

write.csv(
  bind_rows(
    panel_a_data %>%
      transmute(
        panel = "A_baskets",
        basket,
        method,
        estimate,
        q025,
        q975,
        n,
        n_events
      ),
    pooled_band %>%
      transmute(
        panel = "A_pooled_band",
        basket,
        method,
        estimate,
        q025,
        q975,
        n,
        n_events
      ),
    shrinkage_data %>%
      transmute(
        panel = "B_shrinkage",
        basket,
        method = "exnexSurv (Gibbs DA)",
        estimate = estimate_gibbs_exnex,
        q025 = q025_gibbs_exnex,
        q975 = q975_gibbs_exnex,
        n = n_events
      ),
    beta_summary %>%
      transmute(
        panel = "C_age",
        basket = NA_character_,
        method,
        estimate,
        q025,
        q975,
        n = NA_integer_,
        n_events = NA_integer_
      ),
    sigma_summary %>%
      transmute(
        panel = "D_sigma",
        basket = NA_character_,
        method,
        estimate,
        q025,
        q975,
        n = NA_integer_,
        n_events = NA_integer_
      )
  ),
  file.path(export_dir, "FigureData_Article_Fig3_TCGA.csv"),
  row.names = FALSE,
  na = ""
)
cat("\n=== TCGA PanCancer Atlas application ===\n")
cat(
  "Patients:",
  nrow(tcga),
  "| Events:",
  sum(tcga$event),
  "| Censoring:",
  round(1 - mean(tcga$event), 3),
  "\n"
)
print(basket_summary, row.names = FALSE)
cat("\nRuntime (seconds):\n")
print(runtime, row.names = FALSE)
cat("\nPosterior summaries (theta, per method and basket):\n")
print(
  posterior_summary %>%
    filter(parameter == "theta") %>%
    select(method, basket, estimate, q025, q975, rhat, bulk_ess) %>%
    as.data.frame(),
  row.names = FALSE
)

invisible(NULL)
