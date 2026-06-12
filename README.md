# RA-AIPW experiments

This repository contains the final reproducibility code for the RA-AIPW numerical experiments in the manuscript.

The repository is intentionally organized in a **no-folder layout**, so that all files can be uploaded directly through the GitHub web interface.

## Main scripts

| File | Purpose |
|---|---|
| `01_synthetic_weak_overlap_experiment.R` | Reproduces the synthetic weak-overlap experiment and the additional diagnostics reported in the Supplementary Material. |
| `02_synthetic_weak_overlap_diagnostics.R` | Convenience wrapper for checking or regenerating the synthetic diagnostic outputs. |
| `03_nhefs_semisynthetic_experiment.R` | Reproduces the NHEFS-based semi-synthetic experiment in Section 5.2. |
| `00_prepare_nhefs_covariates.R` | Documents the NHEFS data source. The final NHEFS script loads `causaldata::nhefs_complete` directly. |
| `run_all.R` | Runs the synthetic and NHEFS scripts in sequence. |

## Software

The scripts use R and the following packages:

```r
parallel
dplyr
tidyr
ggplot2
readr
tibble
stringr
causaldata
```

Missing packages are installed automatically from CRAN.

## Running the code

For a short smoke test:

```r
Sys.setenv(RA_AIPW_QUICK_TEST = "TRUE")
source("run_all.R")
```

For the full submitted Monte Carlo design:

```r
Sys.unsetenv("RA_AIPW_QUICK_TEST")
source("run_all.R")
```

The full run uses 500 Monte Carlo repetitions in both experiments and may take substantial time.

## Synthetic weak-overlap outputs

The synthetic script produces, among other files,

```text
table_finite_L2grid.csv
table_efficiency_L2grid.csv
fig_rmse_coarse_L2grid.png
fig_tail_q99_coarse_L2grid.png
fig_delta_q95_L2grid.png
fig_L2grid_train_mass.png
adaptive_TS_L2grid_detail.csv
adaptive_TS_L2grid_summary.csv
```

The figures and tables correspond to the additional diagnostics for the synthetic weak-overlap experiment.

## NHEFS-based semi-synthetic outputs

The NHEFS script loads the public `nhefs_complete` data object from the `causaldata` R package and produces

```text
raw_results.csv
raw_results_checkpoint.csv
summary_results.csv
summary_key_results.csv
best_by_kappa.csv
fig_true_propensity_overlap.png
fig_rmse_vs_kappa_color.png
fig_tail_q99_vs_kappa_color.png
fig_coverage_vs_kappa_color.png
fig_active_frac_vs_kappa_color.png
fig_delta_q95_vs_kappa_color.png
```

No individual-level NHEFS data file is distributed in this repository.

## Precomputed outputs

The repository may also include selected CSV and PNG outputs from the submitted runs, so that readers can inspect the numerical results without rerunning the full simulations.

## Code availability statement

```latex
\paragraph{Code availability.}
The code used to reproduce the synthetic weak-overlap diagnostics and the
NHEFS-based semi-synthetic experiment is available at
\url{https://github.com/shinto-eguchi/ra-aipw-experiments}.
```
