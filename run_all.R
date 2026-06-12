############################################################
## Run all final RA-AIPW reproducibility scripts.
##
## Full reproduction uses 500 Monte Carlo repetitions and can take time.
## For a quick smoke test, run first:
##
##   Sys.setenv(RA_AIPW_QUICK_TEST = "TRUE")
##   source("run_all.R")
##
## For the full submitted design, leave RA_AIPW_QUICK_TEST unset or FALSE.
############################################################

source("01_synthetic_weak_overlap_experiment.R")
source("03_nhefs_semisynthetic_experiment.R")

capture.output(sessionInfo(), file = "sessionInfo.txt")
cat("Done. Outputs are written to the repository root.\n")
