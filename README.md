# Reproducibility archive for exnexSurv

This archive contains the analysis scripts, fixed inputs, archived fit summaries, evidence exports, and figures for the article.

> Data-Augmented Gibbs Sampling for EXNEX Log-Normal Survival Models in Basket Trials with Right-Censored Endpoints

The development package lives at <https://github.com/victorney/exnexSurv> and is distributed on CRAN.

## Requirements

The archived runs used R 4.6.0 on Windows 11, rstan 2.32.7 (Stan 2.32.2), and `exnexSurv` 1.3.0 from CRAN, on a 13th-generation Intel Core i7-1355U (10 cores; four used for the chains) with 16 GB of RAM. The scripts also use `survival`, `dplyr`, `tidyr`, `posterior`, `readr`, `tibble`, `ggplot2`, `patchwork`, `purrr`, and `rstan`.

```r
install.packages(c(
  "survival", "dplyr", "tidyr", "posterior", "readr",
  "tibble", "ggplot2", "patchwork", "purrr", "rstan"
))
install.packages("exnexSurv")  # 1.3.0 on CRAN
```

On Windows, `rstan` and `exnexSurv` compile C++ code and may need Rtools.

`exnexSurv` 1.4.0 changed `parallel_chains` from a chain count to a logical flag and deprecated `exnex_surv()` in favour of `pooling_surv()`. Scripts `02`, `05`, `06`, `09`, and `10` select the argument form according to the installed version, so they run against 1.3.x and 1.4.0. The sampler itself is unchanged: refitting one archived Mixed-efficacy replicate under 1.3.0 with long chains agreed with the stored posterior summaries within Monte Carlo error (largest absolute z-score 1.6 across the nine baskets).

## Seed policy

Data generation uses `set.seed(1 + 10000 * seed_index + replicate)`, where `seed_index` is pinned in `01_generate_data.R`: `homogeneous` 1, `mixed` 2, `heterogeneous` 3, `weibull_misspec` 5, `hetvar_misspec` 6, `t_misspec` 7, `infcens_misspec` 8, `tcga_calib` 9. `mixed_divlarge` reuses the `mixed` index, so the two designs share a random-number stream replicate by replicate.

Model fitting uses seeds `1000000 + 10000 * seed_index + 10 * replicate + k`, with `k = 1` for the `exnexSurv` Gibbs sampler and `k = 2, 3, 4` for the Stan EXNEX, Complete-pooling, and No-pooling comparators. The other scripts use their own fixed bases: `05` uses seed 42, `06` uses `3000000 + 10 * replicate + k`, `07` uses `110000 + replicate`, `09` uses `70000 + 10000 * configuration`, and `10` documents its data and fit seeds in the script. Every fit is deterministic from the archived inputs.

## Contents

- `figures/`: `Main_Fig3_Computational.pdf` (article Figure 1), `Main_Fig2_Mixed_Cohorts.pdf` (Figure 2), and `Article_Fig3_TCGA.pdf` (Figure 3).
- `01_generate_data.R`: deterministic data generation for the archived trials.
- `02_fit_scenario.R`: fits one scenario with four chains, 1,000 warm-up and 1,000 retained iterations, for `exnexSurv` and the three Stan comparators. Supports sharding through `FIT_REP_FROM`, `FIT_REP_TO`, and `FIT_SHARD`.
- `03_analyze_article.R`: derives the article evidence exports from the archived fit summaries. It does not fit models.
- `04_generate_article_figures.R`: regenerates the three article figures.
- `05_tcga_application.R`: reproduces the observational TCGA illustration.
- `06_tcga_calibrated_simulation.R`: reproduces the 1,000-trial TCGA-calibrated simulation, with the Stan comparators restricted to the first 100 replicates.
- `07_stan_stress_benchmark.R`: fits the centered Stan model on the 100-trial Complete Heterogeneity stress subset, records divergences and maximum treedepth (export 18), and derives the Gibbs-versus-Stan stress summary (export 22).
- `08_divergent_position.R`: the divergent-position robustness check (Mixed Efficacy versus `mixed_divlarge`, paired on the same replicates, with qualification rates). Requires `02_fit_scenario.R mixed_divlarge` first. Writes exports 26, 27, and 28.
- `09_prior_sensitivity.R`: reproduces `12_prior_sensitivity.csv` by re-fitting the first 500 Mixed-efficacy replicates under nine prior configurations (Gibbs sampler only). Generate the replicates first (see below).
- `10_crossover_benchmark.R`: the basket-size versus censoring crossover behind Table 4. It generates its own datasets and writes exports 29 and 30.
- `models/`: the centered EXNEX Stan model and the Complete-pooling and No-pooling comparators.
- `fitted_models/scenarios/`: archived fit summaries for all scenarios, including the TCGA-calibrated design, the paired divergent-position shards, and the list of Stan subset replicates for the calibrated design.
- `fitted_models/tcga_application/`: fixed TCGA application summaries and posterior-result snapshot. `03_mcmc_diagnostics.csv` reports Rhat and bulk ESS for every parameter of each fit; the accuracy-matched index takes the minimum over the nine `theta` rows, so filter with `parameter == "theta"` before computing it.
- `fitted_models/article_exports/`: the compact exports behind the article tables, figures, and computational diagnostics.
- `data/tcga/raw/`: the nine open TCGA clinical files used in the illustration.

`18_stan_stress_benchmark.csv` and `22_gibbs_on_stan_stress.csv` are regenerated by `07_stan_stress_benchmark.R`. `29_crossover_benchmark.csv` (per dataset) and `30_crossover_benchmark_summary.csv` (index ratios) are regenerated by `10_crossover_benchmark.R`. `31_tcga_calibrated_qualification.csv` is derived by `03_analyze_article.R` from the TCGA-calibrated checkpoint.

The participant-level simulation files (`data/sim_*.rds`) are not archived. The generator recreates them from the seed policy, and the archived checkpoints and compact exports preserve every result used by the article.

## Script outputs

| Script | Outputs |
| --- | --- |
| `01_generate_data.R` | `data/sim_*.rds` (regenerated on demand) |
| `02_fit_scenario.R` | `fitted_models/scenarios/results_scenario_<scenario>[_shard<k>].rds` |
| `03_analyze_article.R` | `fitted_models/article_exports/` exports 01, 02, 04-07, 14, 15, 20, 25, and 31 (no export is numbered 03) |
| `04_generate_article_figures.R` | the three PDFs in `figures/` |
| `05_tcga_application.R` | `fitted_models/tcga_application/` summaries, `article_exports/16_tcga_rmst.csv`, and `FigureData_Article_Fig3_TCGA.csv` |
| `06_tcga_calibrated_simulation.R` | `results_scenario_tcga_calib.rds` and `tcga_calib_stan_subset.csv` |
| `07_stan_stress_benchmark.R` | exports 18 and 22 |
| `08_divergent_position.R` | exports 26, 27, and 28 |
| `09_prior_sensitivity.R` | `12_prior_sensitivity.csv` and `results_prior_sensitivity_*.rds` |
| `10_crossover_benchmark.R` | exports 29 and 30 |

## Reproduction order

From this directory, the quick inspection path uses the supplied checkpoints and does not run any fits:

```text
Rscript 03_analyze_article.R
Rscript 04_generate_article_figures.R
```

The generated CSVs and PDFs match the supplied evidence up to software-dependent CSV formatting and plot rendering details.

Runtime columns (`runtime_seconds` in `29_crossover_benchmark.csv` and the trial-level runtimes behind the `04_*` exports) and the ESS-normalized index `t*` derived from them are wall-clock measurements from the benchmark machine and will differ on other hardware. Deterministic from the seed policy and expected to match exactly are `n_total`, `observed_censoring`, `min_ess`, `mean_ess`, `n_divergent`, and the `29` to `30` aggregation (mean of `t*` per cell, then the Stan/Gibbs ratio).

Regenerating simulation inputs requires `01_generate_data.R`. On Windows PowerShell:

```text
$env:SIM_SCENARIOS = "mixed"; $env:SIM_N_REPS = "500"; Rscript 01_generate_data.R
```

On POSIX shells:

```text
SIM_SCENARIOS=mixed SIM_N_REPS=500 Rscript 01_generate_data.R
```

Then:

- `02_fit_scenario.R <scenario>` fits one scenario. The archived Mixed-efficacy design took about 28 h for 2,000 trials, at roughly 0.72 s per Gibbs fit and 17 s per Stan comparator.
- `07_stan_stress_benchmark.R` needs the first 100 `heterogeneous` replicates and takes about 15 min.
- `09_prior_sensitivity.R` needs the first 500 `mixed` replicates and takes about 1 to 2 h for the nine configurations.
- `10_crossover_benchmark.R` generates its own datasets and takes about 20 min.
- `05_tcga_application.R` uses the included raw files and takes about 2 min.
- `06_tcga_calibrated_simulation.R` fits 1,000 Gibbs replicates and 100 Stan replicates, about 2 h in total.

`02_fit_scenario.R` can be split across processes, which is how the longer scenarios were fitted:

```text
FIT_REP_FROM=1   FIT_REP_TO=120 FIT_SHARD=1 Rscript 02_fit_scenario.R mixed_divlarge
FIT_REP_FROM=121 FIT_REP_TO=240 FIT_SHARD=2 Rscript 02_fit_scenario.R mixed_divlarge
```

Each shard writes `results_scenario_<scenario>_shard<k>.rds` and resumes from it if the process is restarted. `06` and `09` checkpoint every 25 replicates. `05`, `07`, and `10` have no checkpointing and restart from the beginning if interrupted.

The TCGA script uses the immutable cBioPortal datahub commit `0cc9138746c08b304f8dac92c31983e0ef44af1d`. The raw files are included, so it normally does not download them again.

## Citation

See `CITATION.cff`.

## License

See `LICENSE`.
