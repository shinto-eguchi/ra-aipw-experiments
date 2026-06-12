############################################################
## NHEFS semi-synthetic weak-overlap experiment
## Relative L2-grid Ada-AIPW under missing outcomes
## Colab R 7-core one-cell script with progress and ETA
############################################################

## =========================================================
## 0. Basic settings
## =========================================================

set.seed(20260528)

## Set RA_AIPW_QUICK_TEST=TRUE for a short smoke test.
## Default FALSE reproduces the full NHEFS-based Monte Carlo design used for the submitted outputs.
QUICK_TEST <- tolower(Sys.getenv("RA_AIPW_QUICK_TEST", unset = "false")) %in% c("true", "1", "yes", "y")

if (QUICK_TEST) {
  N_REP   <- 20
  K_FOLDS <- 2
  KAPPAS  <- c(0, 2, 3)
} else {
  N_REP   <- 500
  K_FOLDS <- 5
  KAPPAS  <- c(0, 1, 2, 3)
}

N_CORES <- 7
TARGET_RESPONSE <- 0.50
EPS_E  <- 0.005
TRIM_E <- 0.05
FIXED_EC <- 0.05
TAU_MAX <- 0.20
RHO_GRID <- c(0.05, 0.10, 0.20, 0.40, 0.60)
SHRINK_TAU_MIN_COUNT <- 20
RIDGE_M <- 1e-6
RIDGE_B <- 5.0

## Use the repository root for GitHub-friendly output filenames.
OUTDIR <- "."
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("NHEFS Relative L2-grid Ada-AIPW experiment\n")
cat("QUICK_TEST =", QUICK_TEST, "\n")
cat("N_REP      =", N_REP, "\n")
cat("K_FOLDS    =", K_FOLDS, "\n")
cat("N_CORES    =", N_CORES, "\n")
cat("KAPPAS     =", paste(KAPPAS, collapse = ", "), "\n")
cat("RHO_GRID   =", paste(RHO_GRID, collapse = ", "), "\n")
cat("OUTDIR     =", OUTDIR, "\n")
cat("============================================================\n\n")

## =========================================================
## 1. Packages
## =========================================================

pkgs <- c("causaldata", "ggplot2", "parallel", "dplyr", "tidyr", "readr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
library(ggplot2)
library(parallel)
library(dplyr)
library(tidyr)
library(readr)

## =========================================================
## 2. Load NHEFS complete data
## =========================================================

data("nhefs_complete", package = "causaldata")
dat0 <- nhefs_complete

needed <- c(
  "wt82_71",
  "qsmk", "sex", "race", "age", "education",
  "smokeintensity", "smokeyrs",
  "exercise", "active", "wt71"
)
missing_vars <- setdiff(needed, names(dat0))
if (length(missing_vars) > 0) {
  stop("Missing variables in nhefs_complete: ", paste(missing_vars, collapse = ", "))
}

dat <- dat0[, needed]
dat <- dat[complete.cases(dat), ]
dat$Y <- dat$wt82_71

cat_vars <- c("qsmk", "sex", "race", "education", "exercise", "active")
for (v in cat_vars) dat[[v]] <- factor(dat[[v]])

num_vars <- c("age", "smokeintensity", "smokeyrs", "wt71")
for (v in num_vars) dat[[paste0("z_", v)]] <- as.numeric(scale(dat[[v]]))

x_formula <- ~ qsmk + sex + race +
  z_age + I(z_age^2) +
  education +
  z_smokeintensity + I(z_smokeintensity^2) +
  z_smokeyrs + I(z_smokeyrs^2) +
  exercise + active +
  z_wt71 + I(z_wt71^2)

X <- model.matrix(x_formula, data = dat)
if (ncol(X) > 1) {
  keep <- c(TRUE, apply(X[, -1, drop = FALSE], 2, sd) > 1e-10)
  X <- X[, keep, drop = FALSE]
}

Y <- as.numeric(dat$Y)
n <- length(Y)
p <- ncol(X)
psi0 <- mean(Y)

cat("NHEFS complete loaded.\n")
cat("n =", n, " p =", p, "\n")
cat("Full-data target psi0 =", round(psi0, 4), "\n\n")

## =========================================================
## 3. Generate weak-overlap response scores
## =========================================================

X_no_intercept <- X[, -1, drop = FALSE]
X_scaled <- scale(X_no_intercept)

beta_g <- seq(1.0, 0.25, length.out = ncol(X_scaled))
beta_g <- beta_g * rep(c(1, -1), length.out = length(beta_g))
g_raw <- as.vector(X_scaled %*% beta_g)
g_score <- as.numeric(scale(g_raw))

find_intercept <- function(kappa, g, target = 0.5) {
  if (abs(kappa) < 1e-12) return(qlogis(target))
  f <- function(c0) mean(plogis(c0 + kappa * g)) - target
  uniroot(f, interval = c(-30, 30))$root
}

make_e_true <- function(kappa) {
  c0 <- find_intercept(kappa, g_score, TARGET_RESPONSE)
  e <- plogis(c0 + kappa * g_score)
  pmin(pmax(e, 1e-8), 1 - 1e-8)
}

diag_e <- do.call(rbind, lapply(KAPPAS, function(kappa) {
  data.frame(kappa = factor(kappa, levels = KAPPAS), e_true = make_e_true(kappa))
}))
p_overlap <- ggplot(diag_e, aes(x = e_true)) +
  geom_histogram(bins = 40, boundary = 0, fill = "steelblue", color = "white") +
  facet_wrap(~ kappa, nrow = 1) +
  labs(
    title = "True response propensity under increasing weak overlap",
    x = "e_kappa(X)", y = "count"
  ) +
  theme_bw(base_size = 13)
print(p_overlap)
ggsave(file.path(OUTDIR, "fig_true_propensity_overlap.png"),
       p_overlap, width = 11, height = 3.5, dpi = 180)

## =========================================================
## 4. Utility functions
## =========================================================

clip_prob <- function(e, eps = EPS_E) pmin(pmax(e, eps), 1 - eps)

safe_mean <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (length(x) <= 1 || all(is.na(x))) return(NA_real_)
  sd(x, na.rm = TRUE)
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(stats::quantile(x, prob, na.rm = TRUE, names = FALSE, type = 8))
}

ridge_coef <- function(Xtr, ytr, lambda = 1e-6) {
  Xtr <- as.matrix(Xtr)
  ytr <- as.numeric(ytr)
  p <- ncol(Xtr)
  ok <- is.finite(ytr) & apply(Xtr, 1, function(z) all(is.finite(z)))
  Xtr <- Xtr[ok, , drop = FALSE]
  ytr <- ytr[ok]
  if (nrow(Xtr) < 2 || length(unique(ytr)) <= 1 && lambda < 1e-10) {
    b <- rep(0, p)
    b[1] <- mean(ytr, na.rm = TRUE)
    return(b)
  }
  pen <- diag(lambda, p)
  pen[1, 1] <- 0
  out <- tryCatch(
    solve(crossprod(Xtr) + pen, crossprod(Xtr, ytr)),
    error = function(e) {
      qr.solve(crossprod(Xtr) + pen + diag(1e-4, p), crossprod(Xtr, ytr))
    }
  )
  as.numeric(out)
}

safe_glm_binom <- function(Xtr, ytr) {
  ytr <- as.numeric(ytr)
  p <- ncol(Xtr)
  if (length(unique(ytr)) < 2) {
    b <- rep(0, p)
    b[1] <- qlogis(pmin(pmax(mean(ytr), 1e-4), 1 - 1e-4))
    return(b)
  }
  fit <- tryCatch(
    suppressWarnings(glm.fit(x = Xtr, y = ytr, family = binomial(),
                             control = glm.control(maxit = 50))),
    error = function(e) NULL
  )
  if (is.null(fit) || any(!is.finite(fit$coefficients))) {
    b <- rep(0, p)
    b[1] <- qlogis(pmin(pmax(mean(ytr), 1e-4), 1 - 1e-4))
    return(b)
  }
  as.numeric(fit$coefficients)
}

make_folds <- function(n, K) {
  sample(rep(seq_len(K), length.out = n))
}

summarize_phi <- function(phi, method, psi0, active = NA_real_, h_l2 = NA_real_,
                          M_rel = NA_real_, active_train = NA_real_, cutoff = NA_real_,
                          aipw_est = NA_real_) {
  phi <- as.numeric(phi)
  est <- mean(phi)
  se <- sd(phi) / sqrt(length(phi))
  cover <- as.numeric(is.finite(se) && (psi0 >= est - 1.96 * se) && (psi0 <= est + 1.96 * se))
  data.frame(
    method = method,
    est = est,
    err = est - psi0,
    se = se,
    cover = cover,
    q95_abs_phi = safe_quantile(abs(phi), 0.95),
    q99_abs_phi = safe_quantile(abs(phi), 0.99),
    max_abs_phi = max(abs(phi), na.rm = TRUE),
    active_frac = active,
    h_l2 = h_l2,
    M_rel_train = M_rel,
    active_train = active_train,
    cutoff_mean = cutoff,
    delta_sqrt_n_from_aipw = ifelse(is.finite(aipw_est), sqrt(length(phi)) * (est - aipw_est), NA_real_),
    stringsAsFactors = FALSE
  )
}

choose_l2_relative_cutoff <- function(e_train, T_train, tau_max = 0.20, rho = 0.20) {
  e_train <- as.numeric(e_train)
  T_train <- as.numeric(T_train)
  ok <- is.finite(e_train) & is.finite(T_train)
  e_train <- e_train[ok]
  T_train <- T_train[ok]
  if (length(e_train) < 5) return(list(cutoff = 0, M_rel = 0, active = 0))
  denom <- sqrt(mean(T_train^2, na.rm = TRUE))
  if (!is.finite(denom) || denom <= 1e-12) denom <- 1
  cand <- sort(unique(as.numeric(quantile(e_train, probs = seq(0.001, tau_max, length.out = 80),
                                         names = FALSE, type = 8, na.rm = TRUE))))
  cand <- cand[is.finite(cand) & cand > 0]
  if (length(cand) == 0) return(list(cutoff = 0, M_rel = 0, active = 0))
  best <- list(cutoff = 0, M_rel = 0, active = 0)
  for (cc in cand) {
    h <- pmin(1, e_train / cc)
    active <- mean(e_train < cc)
    M_rel <- sqrt(mean(((1 - h)^2) * T_train^2, na.rm = TRUE)) / denom
    if (is.finite(active) && is.finite(M_rel) && active <= tau_max + 1e-12 && M_rel <= rho + 1e-12) {
      if (cc > best$cutoff) best <- list(cutoff = cc, M_rel = M_rel, active = active)
    }
  }
  best
}

## =========================================================
## 5. One Monte Carlo replication
## =========================================================

run_one <- function(job) {
  rep_id <- job$rep_id
  kappa  <- job$kappa
  set.seed(1000003 + 1009 * rep_id + round(100 * kappa))

  e_true <- make_e_true(kappa)
  R <- rbinom(n, size = 1, prob = e_true)
  folds <- make_folds(n, K_FOLDS)

  e_hat <- rep(NA_real_, n)
  m_hat <- rep(NA_real_, n)
  b_hat <- rep(NA_real_, n)

  phi_list <- list(
    OR = rep(NA_real_, n),
    IPW = rep(NA_real_, n),
    AIPW = rep(NA_real_, n),
    TruncAIPW_0.05 = rep(NA_real_, n),
    FixedTS_ec0.05 = rep(NA_real_, n),
    ShrinkQuantileTS = rep(NA_real_, n)
  )
  for (rr in RHO_GRID) phi_list[[paste0("RelAda_rho", sprintf("%.2f", rr))]] <- rep(NA_real_, n)

  diag_by_method <- list()
  for (nm in names(phi_list)) {
    diag_by_method[[nm]] <- list(active = c(), h_l2 = c(), M_rel = c(), active_train = c(), cutoff = c())
  }

  for (k in seq_len(K_FOLDS)) {
    idx_te <- which(folds == k)
    idx_tr <- which(folds != k)
    Xtr <- X[idx_tr, , drop = FALSE]
    Xte <- X[idx_te, , drop = FALSE]
    Rtr <- R[idx_tr]
    Ytr <- Y[idx_tr]

    ## propensity model
    coef_e <- safe_glm_binom(Xtr, Rtr)
    e_tr <- clip_prob(as.vector(plogis(Xtr %*% coef_e)))
    e_te <- clip_prob(as.vector(plogis(Xte %*% coef_e)))

    ## outcome model among observed training cases
    obs_tr <- which(Rtr == 1)
    if (length(obs_tr) >= ncol(Xtr) + 2) {
      coef_m <- ridge_coef(Xtr[obs_tr, , drop = FALSE], Ytr[obs_tr], lambda = RIDGE_M)
    } else if (length(obs_tr) >= 2) {
      coef_m <- ridge_coef(Xtr[obs_tr, , drop = FALSE], Ytr[obs_tr], lambda = 1.0)
    } else {
      coef_m <- rep(0, ncol(Xtr)); coef_m[1] <- mean(Y, na.rm = TRUE)
    }
    m_tr <- as.vector(Xtr %*% coef_m)
    m_te <- as.vector(Xte %*% coef_m)

    ## residual bias-recovery b(X), fitted among observed training cases
    if (length(obs_tr) >= 5) {
      resid_obs <- Ytr[obs_tr] - m_tr[obs_tr]
      coef_b <- ridge_coef(Xtr[obs_tr, , drop = FALSE], resid_obs, lambda = RIDGE_B)
    } else {
      coef_b <- rep(0, ncol(Xtr))
    }
    b_tr <- as.vector(Xtr %*% coef_b)
    b_te <- as.vector(Xte %*% coef_b)

    e_hat[idx_te] <- e_te
    m_hat[idx_te] <- m_te
    b_hat[idx_te] <- b_te

    z_te <- R[idx_te] / e_te * (Y[idx_te] - m_te)
    z_tr <- Rtr / e_tr * (Ytr - m_tr)
    T_tr <- z_tr - b_tr

    ## OR/IPW/AIPW/truncation
    phi_list$OR[idx_te] <- m_te
    phi_list$IPW[idx_te] <- R[idx_te] / e_te * Y[idx_te]
    phi_list$AIPW[idx_te] <- m_te + z_te

    e_trunc_te <- pmax(e_te, TRIM_E)
    phi_list$TruncAIPW_0.05[idx_te] <- m_te + R[idx_te] / e_trunc_te * (Y[idx_te] - m_te)

    ## Fixed TS-AIPW
    h_fixed <- pmin(1, e_te / FIXED_EC)
    phi_list$FixedTS_ec0.05[idx_te] <- m_te + h_fixed * z_te + (1 - h_fixed) * b_te
    diag_by_method$FixedTS_ec0.05$active <- c(diag_by_method$FixedTS_ec0.05$active, mean(h_fixed < 1))
    diag_by_method$FixedTS_ec0.05$h_l2 <- c(diag_by_method$FixedTS_ec0.05$h_l2, sqrt(mean((1 - h_fixed)^2)))

    ## Shrinking quantile TS
    tau_n <- max(n^(-2/3), SHRINK_TAU_MIN_COUNT / n)
    tau_n <- min(tau_n, TAU_MAX)
    c_shrink <- as.numeric(quantile(e_tr, probs = tau_n, names = FALSE, type = 8, na.rm = TRUE))
    if (!is.finite(c_shrink) || c_shrink <= 0) c_shrink <- 0
    h_shrink <- if (c_shrink > 0) pmin(1, e_te / c_shrink) else rep(1, length(idx_te))
    phi_list$ShrinkQuantileTS[idx_te] <- m_te + h_shrink * z_te + (1 - h_shrink) * b_te
    diag_by_method$ShrinkQuantileTS$active <- c(diag_by_method$ShrinkQuantileTS$active, mean(h_shrink < 1))
    diag_by_method$ShrinkQuantileTS$h_l2 <- c(diag_by_method$ShrinkQuantileTS$h_l2, sqrt(mean((1 - h_shrink)^2)))
    diag_by_method$ShrinkQuantileTS$cutoff <- c(diag_by_method$ShrinkQuantileTS$cutoff, c_shrink)

    ## Relative L2-grid Ada-AIPW
    for (rr in RHO_GRID) {
      nm <- paste0("RelAda_rho", sprintf("%.2f", rr))
      ch <- choose_l2_relative_cutoff(e_tr, T_tr, tau_max = TAU_MAX, rho = rr)
      h_ad <- if (ch$cutoff > 0) pmin(1, e_te / ch$cutoff) else rep(1, length(idx_te))
      phi_list[[nm]][idx_te] <- m_te + h_ad * z_te + (1 - h_ad) * b_te
      diag_by_method[[nm]]$active <- c(diag_by_method[[nm]]$active, mean(h_ad < 1))
      diag_by_method[[nm]]$h_l2 <- c(diag_by_method[[nm]]$h_l2, sqrt(mean((1 - h_ad)^2)))
      diag_by_method[[nm]]$M_rel <- c(diag_by_method[[nm]]$M_rel, ch$M_rel)
      diag_by_method[[nm]]$active_train <- c(diag_by_method[[nm]]$active_train, ch$active)
      diag_by_method[[nm]]$cutoff <- c(diag_by_method[[nm]]$cutoff, ch$cutoff)
    }
  }

  aipw_est <- mean(phi_list$AIPW)
  out <- do.call(rbind, lapply(names(phi_list), function(nm) {
    d <- diag_by_method[[nm]]
    summarize_phi(
      phi = phi_list[[nm]], method = nm, psi0 = psi0,
      active = safe_mean(d$active),
      h_l2 = safe_mean(d$h_l2),
      M_rel = safe_mean(d$M_rel),
      active_train = safe_mean(d$active_train),
      cutoff = safe_mean(d$cutoff),
      aipw_est = aipw_est
    )
  }))
  out$rep_id <- rep_id
  out$kappa <- kappa
  out$response_rate <- mean(R)
  out$e_min <- min(e_true)
  out$e_q01 <- safe_quantile(e_true, 0.01)
  out$e_q05 <- safe_quantile(e_true, 0.05)
  out$psi0 <- psi0
  out[, c("rep_id", "kappa", "method", "est", "err", "se", "cover",
          "q95_abs_phi", "q99_abs_phi", "max_abs_phi",
          "active_frac", "h_l2", "M_rel_train", "active_train", "cutoff_mean",
          "delta_sqrt_n_from_aipw", "response_rate", "e_min", "e_q01", "e_q05", "psi0")]
}

## =========================================================
## 6. Run jobs with progress / ETA
##    Colab-safe parallelization by mclapply (forking).
##    This avoids socket-cluster serialization errors such as
##    "invalid connection" during clusterExport/clusterCall.
## =========================================================

jobs <- expand.grid(rep_id = seq_len(N_REP), kappa = KAPPAS, KEEP.OUT.ATTRS = FALSE)
jobs <- jobs[order(jobs$kappa, jobs$rep_id), ]
N_JOBS <- nrow(jobs)
cat("Total jobs =", N_JOBS, "\n")

## Colab is Linux, so parallel::mclapply is available and more stable than PSOCK clusters.
## If forked parallelism is not supported, set N_CORES <- 1 near the top.
USE_MCLAPPLY <- (.Platform$OS.type == "unix" && N_CORES > 1)
if (USE_MCLAPPLY) {
  cat("Parallel backend: mclapply with", min(N_CORES, parallel::detectCores()), "cores\n")
} else {
  cat("Parallel backend: sequential lapply\n")
}

all_raw <- list()
start_time <- Sys.time()
batch_size <- if (USE_MCLAPPLY) min(N_CORES, N_JOBS) else 1
nb <- ceiling(N_JOBS / batch_size)

for (bb in seq_len(nb)) {
  ii <- ((bb - 1) * batch_size + 1):min(bb * batch_size, N_JOBS)
  batch_jobs <- split(jobs[ii, , drop = FALSE], seq_len(length(ii)))
  batch_start <- Sys.time()

  batch_res <- if (USE_MCLAPPLY) {
    parallel::mclapply(
      batch_jobs,
      function(j) {
        ## Each worker gets its own deterministic seed inside run_one().
        run_one(j[1, ])
      },
      mc.cores = min(N_CORES, length(batch_jobs)),
      mc.preschedule = FALSE
    )
  } else {
    lapply(batch_jobs, function(j) run_one(j[1, ]))
  }

  all_raw <- c(all_raw, batch_res)

  done <- length(all_raw)
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  rate <- done / max(elapsed, 1e-8)
  eta <- (N_JOBS - done) / max(rate, 1e-8)
  batch_elapsed <- as.numeric(difftime(Sys.time(), batch_start, units = "mins"))
  cat(sprintf(
    "[progress] %d/%d jobs (%.1f%%) | elapsed %.2f min | ETA %.2f min | last batch %.2f min | %s\n",
    done, N_JOBS, 100 * done / N_JOBS, elapsed, eta, batch_elapsed,
    format(Sys.time(), "%H:%M:%S")
  ))

  if (bb %% 5 == 0 || bb == nb) {
    tmp <- do.call(rbind, all_raw)
    write.csv(tmp, file.path(OUTDIR, "raw_results_checkpoint.csv"), row.names = FALSE)
  }
}

raw <- do.call(rbind, all_raw)
write.csv(raw, file.path(OUTDIR, "raw_results.csv"), row.names = FALSE)

## =========================================================
## 7. Summaries
## =========================================================

summary_results <- raw %>%
  group_by(kappa, method) %>%
  summarise(
    n_rep = n_distinct(rep_id),
    mean_est = mean(est, na.rm = TRUE),
    mean_bias = mean(err, na.rm = TRUE),
    mc_sd = sd(est, na.rm = TRUE),
    mean_se = mean(se, na.rm = TRUE),
    rmse = sqrt(mean(err^2, na.rm = TRUE)),
    coverage = mean(cover, na.rm = TRUE),
    mean_q95_abs_phi = mean(q95_abs_phi, na.rm = TRUE),
    mean_q99_abs_phi = mean(q99_abs_phi, na.rm = TRUE),
    median_max_abs_phi = median(max_abs_phi, na.rm = TRUE),
    active_frac_mean = safe_mean(active_frac),
    h_l2_mean = safe_mean(h_l2),
    M_rel_train = safe_mean(M_rel_train),
    active_train = safe_mean(active_train),
    cutoff_mean = safe_mean(cutoff_mean),
    delta_sd = sd(delta_sqrt_n_from_aipw, na.rm = TRUE),
    delta_mad = median(abs(delta_sqrt_n_from_aipw), na.rm = TRUE),
    delta_q95 = safe_quantile(abs(delta_sqrt_n_from_aipw), 0.95),
    delta_q99 = safe_quantile(abs(delta_sqrt_n_from_aipw), 0.99),
    response_rate = mean(response_rate, na.rm = TRUE),
    e_min = mean(e_min, na.rm = TRUE),
    e_q01 = mean(e_q01, na.rm = TRUE),
    e_q05 = mean(e_q05, na.rm = TRUE),
    .groups = "drop"
  )

method_labels <- c(
  "OR" = "OR",
  "IPW" = "IPW",
  "AIPW" = "AIPW",
  "TruncAIPW_0.05" = "Truncated AIPW",
  "FixedTS_ec0.05" = "Fixed TS",
  "ShrinkQuantileTS" = "Shrinking TS"
)
for (rr in RHO_GRID) method_labels[paste0("RelAda_rho", sprintf("%.2f", rr))] <- paste0("RelAda rho=", sprintf("%.2f", rr))

summary_results$method_label <- unname(method_labels[summary_results$method])
summary_results$method_label <- factor(summary_results$method_label, levels = unname(method_labels))
summary_results <- summary_results %>% arrange(kappa, method_label)

write.csv(summary_results, file.path(OUTDIR, "summary_results.csv"), row.names = FALSE)

key_methods <- c("AIPW", "TruncAIPW_0.05", "FixedTS_ec0.05", "ShrinkQuantileTS",
                 paste0("RelAda_rho", sprintf("%.2f", c(0.10, 0.20, 0.40, 0.60))))
summary_key <- summary_results %>% filter(method %in% key_methods)
write.csv(summary_key, file.path(OUTDIR, "summary_key_results.csv"), row.names = FALSE)

best_by_kappa <- summary_results %>%
  filter(grepl("RelAda", method)) %>%
  group_by(kappa) %>%
  slice_min(order_by = rmse, n = 1, with_ties = FALSE) %>%
  ungroup()
write.csv(best_by_kappa, file.path(OUTDIR, "best_by_kappa.csv"), row.names = FALSE)

cat("\n===== Key summary =====\n")
print(summary_key %>%
        select(kappa, method_label, mean_bias, mc_sd, rmse, coverage,
               mean_q99_abs_phi, median_max_abs_phi, active_frac_mean, M_rel_train, delta_q95),
      n = 200)

cat("\n===== Best Relative Ada by kappa =====\n")
print(best_by_kappa %>% select(kappa, method_label, rmse, coverage, active_frac_mean, M_rel_train), n = 20)

## =========================================================
## 8. Figures
## =========================================================

fig_methods <- c("AIPW", "TruncAIPW_0.05", "FixedTS_ec0.05", "ShrinkQuantileTS",
                 paste0("RelAda_rho", sprintf("%.2f", c(0.10, 0.20, 0.40, 0.60))))
plot_df <- summary_results %>% filter(method %in% fig_methods)

p_rmse <- ggplot(plot_df, aes(x = kappa, y = rmse, color = method_label, group = method_label)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  labs(title = "NHEFS semi-synthetic weak-overlap experiment",
       subtitle = "RMSE for the full-data mean target",
       x = expression(kappa), y = "RMSE", color = "Method") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")
print(p_rmse)
ggsave(file.path(OUTDIR, "fig_rmse_vs_kappa_color.png"), p_rmse, width = 9.5, height = 5.5, dpi = 180)

p_tail <- ggplot(plot_df, aes(x = kappa, y = mean_q99_abs_phi, color = method_label, group = method_label)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  labs(title = "Contribution-tail diagnostic",
       subtitle = "Mean 0.99 quantile of |phi_i| across replications",
       x = expression(kappa), y = expression(q[.99](abs(phi[i]))), color = "Method") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")
print(p_tail)
ggsave(file.path(OUTDIR, "fig_tail_q99_vs_kappa_color.png"), p_tail, width = 9.5, height = 5.5, dpi = 180)

p_cov <- ggplot(plot_df, aes(x = kappa, y = coverage, color = method_label, group = method_label)) +
  geom_hline(yintercept = 0.95, linetype = 2) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  labs(title = "Naive full-target Wald coverage",
       x = expression(kappa), y = "Coverage", color = "Method") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")
print(p_cov)
ggsave(file.path(OUTDIR, "fig_coverage_vs_kappa_color.png"), p_cov, width = 9.5, height = 5.5, dpi = 180)

active_df <- plot_df %>% filter(grepl("TS|RelAda", method))
p_active <- ggplot(active_df, aes(x = kappa, y = active_frac_mean, color = method_label, group = method_label)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  labs(title = "Active stabilization fraction",
       x = expression(kappa), y = "Mean active fraction", color = "Method") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")
print(p_active)
ggsave(file.path(OUTDIR, "fig_active_frac_vs_kappa_color.png"), p_active, width = 9.5, height = 5.5, dpi = 180)

p_delta <- ggplot(plot_df %>% filter(method != "AIPW"),
                  aes(x = kappa, y = delta_q95, color = method_label, group = method_label)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  labs(title = expression(paste("Difference from AIPW: ", q[.95], " of ", abs(sqrt(n)(hat(psi)-hat(psi)[AIPW])))),
       x = expression(kappa), y = expression(q[.95](abs(Delta[n]))), color = "Method") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")
print(p_delta)
ggsave(file.path(OUTDIR, "fig_delta_q95_vs_kappa_color.png"), p_delta, width = 9.5, height = 5.5, dpi = 180)

## =========================================================
## 9. Archive
## =========================================================

writeLines(capture.output(sessionInfo()), con = file.path(OUTDIR, "sessionInfo.txt"))
readme <- c(
  "NHEFS Relative L2-grid Ada-AIPW semi-synthetic weak-overlap experiment",
  "",
  "Main files:",
  "- raw_results.csv",
  "- summary_results.csv",
  "- summary_key_results.csv",
  "- best_by_kappa.csv",
  "- fig_rmse_vs_kappa_color.png",
  "- fig_tail_q99_vs_kappa_color.png",
  "- fig_coverage_vs_kappa_color.png",
  "- fig_active_frac_vs_kappa_color.png",
  "- fig_delta_q95_vs_kappa_color.png",
  "- fig_true_propensity_overlap.png",
  "- sessionInfo.txt"
)
writeLines(readme, con = file.path(OUTDIR, "README.txt"))

## Archive creation is disabled in the GitHub version because all outputs
## are written directly to the repository root.
CREATE_ARCHIVE <- tolower(Sys.getenv("RA_AIPW_CREATE_ARCHIVE", unset = "false")) %in% c("true", "1", "yes", "y")
if (CREATE_ARCHIVE) {
  ZIPFILE <- "nhefs_outputs.zip"
  out_files <- c("raw_results.csv", "raw_results_checkpoint.csv", "summary_results.csv",
                 "summary_key_results.csv", "best_by_kappa.csv",
                 "fig_true_propensity_overlap.png", "fig_rmse_vs_kappa_color.png",
                 "fig_tail_q99_vs_kappa_color.png", "fig_coverage_vs_kappa_color.png",
                 "fig_active_frac_vs_kappa_color.png", "fig_delta_q95_vs_kappa_color.png",
                 "README.txt", "sessionInfo.txt")
  out_files <- out_files[file.exists(out_files)]
  utils::zip(zipfile = ZIPFILE, files = out_files)
}

cat("\n============================================================\n")
cat("Finished NHEFS Relative L2-grid Ada-AIPW experiment\n")
cat("Elapsed total minutes:", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2), "\n")
cat("Output directory:", OUTDIR, "\n")
if (exists("ZIPFILE")) cat("Archive:", ZIPFILE, "\n")
cat("Key CSV:", file.path(OUTDIR, "summary_key_results.csv"), "\n")
cat("============================================================\n")
