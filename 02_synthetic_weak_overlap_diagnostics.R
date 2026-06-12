## Additional diagnostics for the synthetic weak-overlap experiment.
## This script generates Supplementary Figs. S1--S3.

source("ra_aipw_functions.R")


if (!file.exists("synthetic_weak_overlap_raw.csv")) {
  message("Raw synthetic results not found; running scripts/01_synthetic_weak_overlap_experiment.R first.")
  source("01_synthetic_weak_overlap_experiment.R")
}

raw <- read.csv("synthetic_weak_overlap_raw.csv")

diag_df <- raw[grepl("^RA-AIPW", raw$method), ]

## S1: relative deleted-residual mass.
png("supp_synthetic_deleted_mass.png", width = 1800, height = 1200, res = 200)
plot(NA, xlim = range(diag_df$d), ylim = range(diag_df$rel_mass, na.rm = TRUE),
     xlab = "Residual-signal strength d",
     ylab = "Relative deleted-residual mass")
for (meth in unique(diag_df$method)) {
  z <- aggregate(rel_mass ~ d, data = diag_df[diag_df$method == meth, ], mean)
  lines(z$d, z$rel_mass, type = "b", pch = 19)
}
legend("topleft", legend = unique(diag_df$method), lty = 1, pch = 19, bty = "n", cex = 0.8)
dev.off()

## S2: upper tail of the stabilized pseudo-outcome.
png("supp_synthetic_tail_q99.png", width = 1800, height = 1200, res = 200)
plot(NA, xlim = range(raw$d), ylim = range(raw$q99_abs_phi, na.rm = TRUE),
     xlab = "Residual-signal strength d",
     ylab = expression(q[0.99](abs(phi[i]))))
for (meth in unique(raw$method)) {
  z <- aggregate(q99_abs_phi ~ d, data = raw[raw$method == meth, ], mean)
  lines(z$d, z$q99_abs_phi, type = "b", pch = 19)
}
legend("topleft", legend = unique(raw$method), lty = 1, pch = 19, bty = "n", cex = 0.7)
dev.off()

## S3: distance from the ordinary AIPW endpoint.
aipw <- raw[raw$method == "AIPW", c("d", "rep", "estimate")]
names(aipw)[3] <- "estimate_aipw"

merged <- merge(raw, aipw, by = c("d", "rep"), all.x = TRUE)
merged$distance_from_aipw <- abs(merged$estimate - merged$estimate_aipw)

png("supp_synthetic_aipw_distance.png", width = 1800, height = 1200, res = 200)
plot(NA, xlim = range(merged$d), ylim = range(merged$distance_from_aipw, na.rm = TRUE),
     xlab = "Residual-signal strength d",
     ylab = "Mean absolute distance from AIPW")
for (meth in unique(merged$method)) {
  if (meth == "AIPW") next
  z <- aggregate(distance_from_aipw ~ d, data = merged[merged$method == meth, ], mean)
  lines(z$d, z$distance_from_aipw, type = "b", pch = 19)
}
legend("topleft", legend = setdiff(unique(merged$method), "AIPW"),
       lty = 1, pch = 19, bty = "n", cex = 0.7)
dev.off()

message("Wrote Supplementary diagnostic figures to .")
