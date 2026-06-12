## Run all reproducibility scripts for the RA-AIPW experiments.


source("ra_aipw_functions.R")

message("Running synthetic weak-overlap experiment...")
source("01_synthetic_weak_overlap_experiment.R")

message("Running additional synthetic diagnostics...")
source("02_synthetic_weak_overlap_diagnostics.R")

message("Preparing NHEFS covariates, if available...")
source("00_prepare_nhefs_covariates.R")

message("Running NHEFS-based semi-synthetic experiment...")
source("03_nhefs_semisynthetic_experiment.R")

capture.output(sessionInfo(), file = "sessionInfo.txt")

message("Done. Results are in  and .")
