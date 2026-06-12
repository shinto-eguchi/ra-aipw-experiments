# RA-AIPW experiments

This is the **no-folder layout** for direct GitHub upload. All files are placed in the repository root.

This repository contains reproducibility materials for the numerical experiments in the paper

**Residual-allocated augmented inverse-probability weighting under weak overlap**

The code is organized around two parts of the manuscript:

1. the synthetic weak-overlap experiment, including additional diagnostics; and  
2. the NHEFS-based semi-synthetic experiment.

## Repository structure

```text
README.md
run_all.R
R/
  ra_aipw_functions.R
scripts/
  00_prepare_nhefs_covariates.R
  01_synthetic_weak_overlap_experiment.R
  02_synthetic_weak_overlap_diagnostics.R
  03_nhefs_semisynthetic_experiment.R


data/
  README.md
sessionInfo.txt
```

## Reproducing the results

From the top-level directory, run

```r
source("run_all.R")
```

This creates the `` and `` directories if needed and runs the synthetic and NHEFS-based scripts in order.

The synthetic scripts are self-contained. The NHEFS-based script expects a prepared covariate file; see `DATA_README.md` and `scripts/00_prepare_nhefs_covariates.R`.

## Manuscript outputs

The scripts are intended to reproduce the following outputs.

| Manuscript item | Script | Output |
|---|---|---|
| Synthetic weak-overlap experiment | `scripts/01_synthetic_weak_overlap_experiment.R` | `synthetic_weak_overlap_summary.csv` |
| Additional diagnostics, Supplementary Figs. S1--S3 | `scripts/02_synthetic_weak_overlap_diagnostics.R` | `supp_synthetic_deleted_mass.png`, `supp_synthetic_tail_q99.png`, `supp_synthetic_aipw_distance.png` |
| NHEFS-based semi-synthetic experiment, Section 5.2 | `scripts/03_nhefs_semisynthetic_experiment.R` | `nhefs_semisynthetic_summary.csv`, `fig_rmse_vs_kappa_color.png`, `fig_tail_q99_vs_kappa_color.png` |

## Data

The synthetic weak-overlap experiment does not require external data.

The NHEFS-based semi-synthetic experiment uses the covariate distribution from the NHANES I Epidemiologic Followup Study (NHEFS). No individual-level data are distributed in this repository. Place the prepared covariate file at

```text
nhefs_covariates.csv
```

or set the environment variable

```text
NHEFS_CSV=/path/to/nhefs_covariates.csv
```

The file should contain one row per subject and columns containing the covariates used to construct the semi-synthetic response and outcome mechanisms. Numeric and factor covariates are handled by `scripts/00_prepare_nhefs_covariates.R`.

## Software

The scripts use base R and recommended packages only (`stats`, `utils`, `graphics`, and `splines`). Package and R-version information are saved to `sessionInfo.txt` when `run_all.R` is executed.

## Random seeds

All simulation scripts use fixed random seeds. The defaults are chosen for reproducibility rather than speed. For quick checks, reduce `n_reps` inside the corresponding script.

## Notes for revision

The current files provide a clean reproducibility scaffold. If the exact Monte Carlo code used to produce the submitted numerical values is available separately, replace the model-specific blocks in `scripts/01_synthetic_weak_overlap_experiment.R` and `scripts/03_nhefs_semisynthetic_experiment.R` while keeping the same output filenames.
