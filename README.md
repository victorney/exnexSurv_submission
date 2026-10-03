# Reproducibility archive for exnexSurv

Analysis scripts, fixed inputs, archived fit summaries, evidence exports, and figures for the article:

> Data-Augmented Gibbs Sampling for EXNEX Log-Normal Survival Models in Basket Trials with Right-Censored Endpoints

The package is on CRAN (`exnexSurv` 1.3.0); development sources are at <https://github.com/victorney/exnexSurv>.

## Requirements

The archived runs used R 4.6.0 on Windows 11, rstan 2.32.7 (Stan 2.32.2), and `exnexSurv` 1.3.0, on a 13th-generation Intel Core i7-1355U (10 cores; four used for the chains) with 16 GB of RAM. The scripts also use `survival`, `dplyr`, `tidyr`, `posterior`, `readr`, `tibble`, `ggplot2`, `patchwork`, `purrr`, and `rstan`.

```r
install.packages(c(
  "survival", "dplyr", "tidyr", "posterior", "readr",
  "tibble", "ggplot2", "patchwork", "purrr", "rstan"
))
install.packages("exnexSurv")  # 1.3.0 on CRAN
```

On Windows, `rstan` and `exnexSurv` compile C++ code and may need Rtools.

## Seed policy

Data generation uses `set.seed(1 + 10000 * seed_index + replicate)`, with `seed_index` pinned in `01_generate_data.R`: `homogeneous` 1, `mixed` 2, `heterogeneous` 3, `weibull_misspec` 5, `hetvar_misspec` 6, `t_misspec` 7, `infcens_misspec` 8, `tcga_calib` 9. `mixed_divlarge` reuses the `mixed` index.

Model fitting uses seeds `1000000 + 10000 * seed_index + 10 * replicate + k`, with `k = 1` for the `exnexSurv` Gibbs sampler and `k = 2, 3, 4` for the Stan comparators. Scripts `05`, `06`, `07`, and `09` use their own fixed bases (seed 42; `3000000 + 10 * replicate + k`; `110000 + replicate`; `70000 + 10000 * configuration`), and `10` documents its seeds in the script. Every fit is deterministic from the archived inputs.

## Contents

- `figures/`: the three article figures (Figure 1: `Main_Fig3_Computational.pdf`; Figure 2: `Main_Fig2_Mixed_Cohorts.pdf`; Figure 3: `Article_Fig3_TCGA.pdf`).
- `01_generate_data.R` through `10_crossover_benchmark.R`: the pipeline, in order (see Script outputs).
- `models/`: the centered EXNEX Stan model and the Complete-pooling and No-pooling comparators.
- `fitted_models/scenarios/`: archived fit summaries, including the TCGA-calibrated design, the divergent-position shards, and the Stan subset list.
- `fitted_models/tcga_application/`: TCGA summaries and posterior snapshot.
- `fitted_models/article_exports/`: the compact exports behind the article tables and figures.
- `data/tcga/raw/`: the nine TCGA clinical files.

`data/sim_*.rds` files are not archived; `01_generate_data.R` recreates them from the seed policy.

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

Quick inspection (no model fitting):

```text
Rscript 03_analyze_article.R
Rscript 04_generate_article_figures.R
```

CSV and PDF output matches the archived evidence up to CSV formatting and plot rendering details. Runtime columns and the `t*` index derived from them are wall-clock measurements from the benchmark machine and will differ on other hardware; `n_total`, censoring, ESS, and the `29`-to-`30` aggregation match exactly.

Regenerating inputs and refitting takes longer: about 28 h for the 2,000-trial Mixed-efficacy design (roughly 0.72 s per Gibbs fit, 17 s per Stan fit), about 2 h for the TCGA-calibrated simulation, 20 min for the crossover, 2 min for the TCGA application. `07` needs the first 100 `heterogeneous` replicates; `09` needs the first 500 `mixed` replicates. `02` supports sharding (`FIT_REP_FROM`, `FIT_REP_TO`, `FIT_SHARD`); `06` and `09` checkpoint every 25 replicates.

Generate inputs on Windows PowerShell:

```text
$env:SIM_SCENARIOS = "mixed"; $env:SIM_N_REPS = "500"; Rscript 01_generate_data.R
```

On POSIX shells:

```text
SIM_SCENARIOS=mixed SIM_N_REPS=500 Rscript 01_generate_data.R
```

The TCGA script uses cBioPortal datahub commit `0cc9138746c08b304f8dac92c31983e0ef44af1d`; the raw files are included, so nothing is downloaded again.

## Citation

See `CITATION.cff`.

## License

See `LICENSE`.
