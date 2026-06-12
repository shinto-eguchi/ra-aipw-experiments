## Common functions for RA-AIPW experiments.
## The implementation uses base R only.

expit <- function(x) {
  1 / (1 + exp(-x))
}

logit <- function(p) {
  log(p / (1 - p))
}

safe_divide <- function(num, den, eps = 1e-8) {
  num / pmax(den, eps)
}

calibrate_intercept <- function(kappa, target = 0.5, x_sd = 1) {
  f <- function(a) {
    mean(expit(a - kappa * rnorm(200000, sd = x_sd))) - target
  }
  ## Use a deterministic quadrature-like Monte Carlo grid for stability.
  z <- qnorm(seq(0.0005, 0.9995, length.out = 20000), sd = x_sd)
  g <- function(a) mean(expit(a - kappa * z)) - target
  uniroot(g, lower = -20, upper = 20)$root
}

make_folds <- function(n, k = 5, seed = 1) {
  set.seed(seed)
  sample(rep(seq_len(k), length.out = n))
}

cap_weight <- function(ehat, cutoff) {
  pmin(1, safe_divide(ehat, cutoff))
}

hard_weight <- function(ehat, cutoff) {
  as.numeric(ehat >= cutoff)
}

observed_residual_contribution <- function(y, r, ehat, mhat) {
  safe_divide(r * (y - mhat), ehat)
}

aipw_pseudo <- function(y, r, ehat, mhat) {
  mhat + observed_residual_contribution(y, r, ehat, mhat)
}

weighted_pseudo <- function(y, r, ehat, mhat, h, completion = 0) {
  mhat + h * observed_residual_contribution(y, r, ehat, mhat) +
    (1 - h) * completion
}

estimate_from_pseudo <- function(phi, psi_true = NA_real_) {
  psi_hat <- mean(phi)
  se <- sd(phi) / sqrt(length(phi))
  out <- list(
    estimate = psi_hat,
    se = se,
    q99_abs_phi = as.numeric(quantile(abs(phi), probs = 0.99, names = FALSE))
  )
  if (!is.na(psi_true)) {
    out$bias <- psi_hat - psi_true
    out$covered <- as.numeric(abs(psi_hat - psi_true) <= 1.96 * se)
  }
  out
}

summarise_mc <- function(df) {
  aggregate_values <- function(x) {
    c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
  }

  methods <- unique(df$method)
  ans <- lapply(methods, function(meth) {
    z <- df[df$method == meth, , drop = FALSE]
    data.frame(
      method = meth,
      bias = mean(z$estimate - z$psi_true, na.rm = TRUE),
      rmse = sqrt(mean((z$estimate - z$psi_true)^2, na.rm = TRUE)),
      coverage = mean(z$covered, na.rm = TRUE),
      q99_abs_phi = mean(z$q99_abs_phi, na.rm = TRUE),
      active = mean(z$active, na.rm = TRUE),
      rel_mass = mean(z$rel_mass, na.rm = TRUE)
    )
  })
  do.call(rbind, ans)
}

choose_cap_by_tolerance <- function(ehat, deleted_basis, tau = 0.20, rho = 0.20,
                                    probs = seq(0.01, 0.60, by = 0.01)) {
  ## Candidate cutoffs are propensity quantiles.
  cand <- unique(as.numeric(quantile(ehat, probs = probs, names = FALSE)))
  cand <- cand[is.finite(cand) & cand > 0]

  total_mass <- mean(abs(deleted_basis)) + 1e-12

  tab <- lapply(cand, function(cutoff) {
    h <- cap_weight(ehat, cutoff)
    active <- mean(h < 1)
    rel_mass <- mean(abs((1 - h) * deleted_basis)) / total_mass
    data.frame(cutoff = cutoff, active = active, rel_mass = rel_mass)
  })
  tab <- do.call(rbind, tab)

  ok <- tab$active <= tau & tab$rel_mass <= rho
  if (!any(ok)) {
    j <- which.min(tab$active + tab$rel_mass)
  } else {
    ## Choose the strongest stabilization among rules satisfying the tolerances.
    j <- max(which(ok))
  }

  list(cutoff = tab$cutoff[j], diagnostics = tab[j, , drop = FALSE], grid = tab)
}

fit_synthetic_nuisance <- function(x, y, r, folds, outcome_formula = y ~ x) {
  n <- length(x)
  ehat <- rep(NA_real_, n)
  mhat <- rep(NA_real_, n)
  completion <- rep(NA_real_, n)

  for (fold in sort(unique(folds))) {
    train <- folds != fold
    test <- folds == fold

    efit <- glm(r ~ x, family = binomial(), subset = train)
    ehat[test] <- predict(efit, newdata = data.frame(x = x[test]), type = "response")

    ofit <- lm(outcome_formula, data = data.frame(y = y, x = x), subset = train & r == 1)
    mhat[test] <- predict(ofit, newdata = data.frame(x = x[test]))

    train_resid <- y[train] - predict(ofit, newdata = data.frame(x = x[train]))
    comp_data <- data.frame(
      resid = train_resid,
      x = x[train],
      xpos = pmax(x[train], 0),
      xpos_ind = as.numeric(x[train] > 0)
    )
    cfit <- lm(resid ~ xpos_ind + xpos, data = comp_data, subset = r[train] == 1)
    completion[test] <- predict(cfit, newdata = data.frame(
      x = x[test],
      xpos = pmax(x[test], 0),
      xpos_ind = as.numeric(x[test] > 0)
    ))
  }

  list(ehat = pmin(pmax(ehat, 1e-4), 1 - 1e-4),
       mhat = mhat,
       completion = completion)
}

simulate_synthetic_data <- function(n = 12000, kappa = 3.5, d = 1,
                                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  x <- rnorm(n)
  a <- calibrate_intercept(kappa, target = 0.5)
  e <- expit(a - kappa * x)
  r <- rbinom(n, size = 1, prob = e)

  mu <- d * (0.35 * as.numeric(x > 0) + 0.75 * pmax(x, 0))
  y_full <- mu + rnorm(n)
  y_obs <- ifelse(r == 1, y_full, 0)

  list(
    x = x,
    r = r,
    y = y_obs,
    y_full = y_full,
    e = e,
    mu = mu,
    psi_true = mean(y_full)
  )
}

evaluate_methods <- function(dat, tau = 0.20, rho_values = c(0.10, 0.20, 0.40),
                             cutoff_fixed = 0.05, folds = NULL) {
  n <- length(dat$x)
  if (is.null(folds)) folds <- make_folds(n, k = 5, seed = 1)

  nuis <- fit_synthetic_nuisance(dat$x, dat$y, dat$r, folds)
  ehat <- nuis$ehat
  mhat <- nuis$mhat
  comp <- nuis$completion

  out <- list()

  phi_aipw <- aipw_pseudo(dat$y, dat$r, ehat, mhat)
  est <- estimate_from_pseudo(phi_aipw, dat$psi_true)
  out[[length(out) + 1]] <- data.frame(
    method = "AIPW",
    estimate = est$estimate,
    se = est$se,
    q99_abs_phi = est$q99_abs_phi,
    covered = est$covered,
    psi_true = dat$psi_true,
    active = NA_real_,
    rel_mass = NA_real_
  )

  h_tr <- hard_weight(ehat, cutoff_fixed)
  phi_tr <- weighted_pseudo(dat$y, dat$r, ehat, mhat, h_tr, completion = 0)
  est <- estimate_from_pseudo(phi_tr, dat$psi_true)
  out[[length(out) + 1]] <- data.frame(
    method = "Truncated AIPW",
    estimate = est$estimate,
    se = est$se,
    q99_abs_phi = est$q99_abs_phi,
    covered = est$covered,
    psi_true = dat$psi_true,
    active = mean(h_tr < 1),
    rel_mass = NA_real_
  )

  h_cap <- cap_weight(ehat, cutoff_fixed)
  phi_cap <- weighted_pseudo(dat$y, dat$r, ehat, mhat, h_cap, completion = 0)
  est <- estimate_from_pseudo(phi_cap, dat$psi_true)
  out[[length(out) + 1]] <- data.frame(
    method = "Fixed cap",
    estimate = est$estimate,
    se = est$se,
    q99_abs_phi = est$q99_abs_phi,
    covered = est$covered,
    psi_true = dat$psi_true,
    active = mean(h_cap < 1),
    rel_mass = NA_real_
  )

  deleted_basis <- observed_residual_contribution(dat$y, dat$r, ehat, mhat)
  for (rho in rho_values) {
    sel <- choose_cap_by_tolerance(ehat, deleted_basis, tau = tau, rho = rho)
    h <- cap_weight(ehat, sel$cutoff)
    phi <- weighted_pseudo(dat$y, dat$r, ehat, mhat, h, completion = comp)
    est <- estimate_from_pseudo(phi, dat$psi_true)
    out[[length(out) + 1]] <- data.frame(
      method = sprintf("RA-AIPW rho=%.2f", rho),
      estimate = est$estimate,
      se = est$se,
      q99_abs_phi = est$q99_abs_phi,
      covered = est$covered,
      psi_true = dat$psi_true,
      active = sel$diagnostics$active,
      rel_mass = sel$diagnostics$rel_mass
    )
  }

  do.call(rbind, out)
}

plot_metric <- function(summary_df, xvar, yvar, groupvar, file, ylab, xlab = xvar) {
  png(file, width = 1800, height = 1200, res = 200)
  on.exit(dev.off(), add = TRUE)

  groups <- unique(summary_df[[groupvar]])
  xs <- sort(unique(summary_df[[xvar]]))
  yr <- range(summary_df[[yvar]], na.rm = TRUE)

  plot(xs, rep(NA_real_, length(xs)), ylim = yr, xlab = xlab, ylab = ylab,
       type = "n")
  for (g in groups) {
    z <- summary_df[summary_df[[groupvar]] == g, , drop = FALSE]
    z <- z[order(z[[xvar]]), , drop = FALSE]
    lines(z[[xvar]], z[[yvar]], type = "b", pch = 19)
  }
  legend("topleft", legend = groups, lty = 1, pch = 19, bty = "n", cex = 0.8)
}
