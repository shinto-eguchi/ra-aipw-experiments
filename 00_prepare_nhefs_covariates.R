############################################################
## NHEFS data note
##
## The final NHEFS semi-synthetic code uses the public data object
## causaldata::nhefs_complete directly.  Therefore no individual-level
## NHEFS file is distributed with this repository and no manual data
## preparation is required.
##
## To run the NHEFS experiment:
##   source("03_nhefs_semisynthetic_experiment.R")
############################################################

if (!requireNamespace("causaldata", quietly = TRUE)) {
  install.packages("causaldata", repos = "https://cloud.r-project.org")
}

data("nhefs_complete", package = "causaldata")
cat("Loaded causaldata::nhefs_complete with", nrow(nhefs_complete), "rows.\n")
