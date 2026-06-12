## Synthetic weak-overlap experiment.
## This script generates synthetic_weak_overlap_summary.csv.

source("ra_aipw_functions.R")


set.seed(20260613)

n <- 12000
kappa <- 3.5
d_grid <- c(0.5, 1, 2, 4)
rho_values <- c(0.10, 0.20, 0.40)
tau <- 0.20
n_reps <- 100

all_results <- list()
counter <- 1

for (d in d_grid) {
  message("Synthetic experiment: d = ", d)

  for (rep in seq_len(n_reps)) {
    dat <- simulate_synthetic_data(
      n = n,
      kappa = kappa,
      d = d,
      seed = 100000 + 1000 * match(d, d_grid) + rep
    )
    folds <- make_folds(n, k = 5, seed = 200000 + rep)

    z <- evaluate_methods(
      dat,
      tau = tau,
      rho_values = rho_values,
      cutoff_fixed = 0.05,
      folds = folds
    )
    z$d <- d
    z$rep <- rep
    all_results[[counter]] <- z
    counter <- counter + 1
  }
}

mc <- do.call(rbind, all_results)

summary_list <- lapply(split(mc, list(mc$d, mc$method), drop = TRUE), function(z) {
  data.frame(
    d = unique(z$d),
    method = unique(z$method),
    bias = mean(z$estimate - z$psi_true, na.rm = TRUE),
    rmse = sqrt(mean((z$estimate - z$psi_true)^2, na.rm = TRUE)),
    coverage = mean(z$covered, na.rm = TRUE),
    q99_abs_phi = mean(z$q99_abs_phi, na.rm = TRUE),
    active = mean(z$active, na.rm = TRUE),
    rel_mass = mean(z$rel_mass, na.rm = TRUE)
  )
})

summary_df <- do.call(rbind, summary_list)
summary_df <- summary_df[order(summary_df$d, summary_df$method), ]

write.csv(mc, "synthetic_weak_overlap_raw.csv", row.names = FALSE)
write.csv(summary_df, "synthetic_weak_overlap_summary.csv", row.names = FALSE)

message("Wrote synthetic_weak_overlap_summary.csv")
