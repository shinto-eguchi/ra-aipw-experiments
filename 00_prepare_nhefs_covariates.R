## Prepare the NHEFS covariate matrix.
##
## This repository does not distribute individual-level NHEFS data.
## Place a prepared CSV file at nhefs_covariates.csv, or set
## the environment variable NHEFS_CSV to its path.
##
## The expected object is a data frame with one row per subject and covariate
## columns only. Numeric columns are standardized; non-numeric columns are
## converted into dummy variables by model.matrix.


input_path <- Sys.getenv("NHEFS_CSV", unset = "nhefs_covariates.csv")
output_path <- "nhefs_covariates_prepared.csv"

prepare_covariate_matrix <- function(dat) {
  dat <- as.data.frame(dat)
  dat <- dat[complete.cases(dat), , drop = FALSE]

  if (nrow(dat) == 0) {
    stop("No complete rows remain after removing missing values.")
  }

  mm <- model.matrix(~ . , data = dat)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]

  is_num <- apply(mm, 2, is.numeric)
  mm[, is_num] <- scale(mm[, is_num, drop = FALSE])

  as.data.frame(mm)
}

if (file.exists(input_path)) {
  message("Reading NHEFS covariates from ", input_path)
  dat <- read.csv(input_path)
  xdat <- prepare_covariate_matrix(dat)
  write.csv(xdat, output_path, row.names = FALSE)
  message("Wrote ", output_path)
} else {
  message("NHEFS covariate file not found.")
  message("Expected: ", input_path)
  message("Skipping preparation. See DATA_README.md.")
}
