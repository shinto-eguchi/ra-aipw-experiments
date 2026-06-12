############################################################
## Adaptive / L2-grid TS-AIPW simulation
##
## Colab R one-cell script, designed for 7 cores.
##
## New finite-sample tuning idea:
##   Instead of using a shrinking quantile cutoff only, choose the
##   cutoff c by the fold-wise training constraints
##
##     tau_hat(c) <= tau_max,
##     M_hat(c)   <= rho,
##
##   where
##
##     h_c(x) = min{1, (e_hat(x)/c)^a},
##     T_hat  = R/e_hat(X) {Y - m_hat(X)} - r_hat(X),
##     M_hat(c)^2 = mean_train[{1-h_c(X)}^2 T_hat^2].
##
##   We try tau_max = 0.20 and rho in {0.05, 0.10, 0.20, 0.40}.
##   This lets the L2 constraint, rather than a very small shrinking
##   quantile, decide how much stabilization is safe.
############################################################

## ==========================================================
## 0. Settings
## ==========================================================

set.seed(20260527)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

## Set RA_AIPW_QUICK_TEST=TRUE for a short smoke test.
## Default FALSE reproduces the full Monte Carlo design used for the submitted outputs.
QUICK_TEST <- tolower(Sys.getenv("RA_AIPW_QUICK_TEST", unset = "false")) %in% c("true", "1", "yes", "y")

N_REP   <- if (QUICK_TEST) 40 else 500
N_CORES <- 7

NSAMP <- if (QUICK_TEST) c(500, 1500) else c(500, 1500, 5000, 12000)
KAPPAS <- c(3.5)
D_EFFECTS <- if (QUICK_TEST) c(1, 4) else c(0.5, 1, 2, 4)
SCENARIOS <- c("correct_m", "coarse_m")

K_FOLD <- 5
SIGMA <- 1
TARGET_RESPONSE <- 0.50

## Benchmarks.
TRIM_E       <- 0.05
FIXED_EC_TS  <- 0.05

## Shrinking quantile benchmark, retained only as an asymptotic-recovery diagnostic.
TAU_POWER  <- 2/3
MIN_ACTIVE <- 10

## New L2-grid tuning.
TAU_MAX  <- 0.20
RHO_GRID <- c(0.05, 0.10, 0.20, 0.40)

## Stabilizer powers.
## a = 1 gives h/e <= 1/c, an effective inverse-weight cap.
## q4 = 3/4 is included as an optional smoother moment-tail rule.
## For a quicker main run, set A_POWERS <- c(a1 = 1).
A_POWERS <- c(a1 = 1)
# A_POWERS <- c(a1 = 1, q4 = 3/4)

EPS_E <- 1e-6
## Use the repository root for GitHub-friendly output filenames.
OUTDIR <- "."
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("Adaptive / L2-grid TS-AIPW simulation\n")
cat("N_REP      =", N_REP, "\n")
cat("N_CORES    =", N_CORES, "\n")
cat("NSAMP      =", paste(NSAMP, collapse = ", "), "\n")
cat("KAPPAS     =", paste(KAPPAS, collapse = ", "), "\n")
cat("D_EFFECTS  =", paste(D_EFFECTS, collapse = ", "), "\n")
cat("SCENARIOS  =", paste(SCENARIOS, collapse = ", "), "\n")
cat("K_FOLD     =", K_FOLD, "\n")
cat("TAU_MAX    =", TAU_MAX, "\n")
cat("RHO_GRID   =", paste(RHO_GRID, collapse = ", "), "\n")
cat("A_POWERS   =", paste(names(A_POWERS), A_POWERS, sep = "=", collapse = ", "), "\n")
cat("OUTDIR     =", OUTDIR, "\n")
cat("============================================================\n\n")

## ==========================================================
## 1. Packages
## ==========================================================

pkgs <- c("parallel", "dplyr", "tidyr", "ggplot2", "readr", "tibble", "stringr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(tibble)
library(stringr)

## ==========================================================
## 2. DGP and utility functions
## ==========================================================

safe_plogis <- function(x) pmin(pmax(plogis(x), EPS_E), 1 - EPS_E)
clamp_e <- function(e) pmin(pmax(as.numeric(e), EPS_E), 1 - EPS_E)

normal_expectation_plogis <- function(a0, kappa) {
  f <- function(x) plogis(a0 - kappa * x) * dnorm(x)
  integrate(f, lower = -Inf, upper = Inf, rel.tol = 1e-10)$value
}

find_a0 <- function(kappa, target = 0.5) {
  f <- function(a0) normal_expectation_plogis(a0, kappa) - target
  uniroot(f, interval = c(-50, 50))$root
}

A0_TABLE <- tibble(
  kappa = KAPPAS,
  a0 = vapply(KAPPAS, find_a0, numeric(1), target = TARGET_RESPONSE)
)

get_a0 <- function(kappa) A0_TABLE$a0[match(kappa, A0_TABLE$kappa)]

e_fun <- function(x, kappa) safe_plogis(get_a0(kappa) - kappa * x)
w_fun <- function(x) as.numeric(x > 0)

mu_fun <- function(x, d_effect) {
  d_effect * (0.35 * w_fun(x) + 0.75 * pmax(x, 0))
}

## Analytic target E{mu(X)} for X ~ N(0,1).
psi_true <- function(d_effect) {
  d_effect * (0.35 * 0.5 + 0.75 / sqrt(2 * pi))
}

tau_shrink_fun <- function(n) max(n^(-TAU_POWER), MIN_ACTIVE / n)

h_stabilizer <- function(e, ec, a = 1) {
  ec <- pmax(as.numeric(ec), EPS_E)
  h <- (pmax(as.numeric(e), EPS_E) / ec)^a
  pmin(1, pmax(0, h))
}

basis_m <- function(x, scenario) {
  if (scenario == "correct_m") {
    ## Contains the true mean form: 1, W, X_+.
    cbind(Intercept = rep(1, length(x)), W = w_fun(x), Xplus = pmax(x, 0))
  } else if (scenario == "coarse_m") {
    ## Deliberately misspecified outcome regression.
    matrix(1, nrow = length(x), ncol = 1, dimnames = list(NULL, "Intercept"))
  } else {
    stop("Unknown scenario: ", scenario)
  }
}

basis_r <- function(x) {
  ## Coarse residual-centering basis. Captures the step part but not X_+.
  cbind(Intercept = rep(1, length(x)), W = w_fun(x))
}

ridge_coef <- function(B, y, ridge = 1e-8) {
  B <- as.matrix(B)
  y <- as.numeric(y)
  row_ok <- is.finite(y) & apply(B, 1, function(z) all(is.finite(z)))
  if (sum(row_ok) == 0) return(rep(0, ncol(B)))
  B0 <- B[row_ok, , drop = FALSE]
  y0 <- y[row_ok]
  if (nrow(B0) < ncol(B0)) {
    ans <- rep(0, ncol(B0))
    ans[1] <- mean(y0, na.rm = TRUE)
    return(ans)
  }
  A <- crossprod(B0) + ridge * diag(ncol(B0))
  b <- crossprod(B0, y0)
  as.numeric(tryCatch(solve(A, b), error = function(e) qr.solve(A, b)))
}

fit_predict_e <- function(x_train, r_train, x_test, x_cutoff_train) {
  dat <- data.frame(R = r_train, X = x_train)
  fit <- tryCatch(
    suppressWarnings(glm(R ~ X, data = dat, family = binomial())),
    error = function(e) NULL
  )
  if (is.null(fit) || any(!is.finite(coef(fit)))) {
    p0 <- mean(r_train)
    e_test  <- rep(p0, length(x_test))
    e_train <- rep(p0, length(x_cutoff_train))
  } else {
    e_test <- tryCatch(
      predict(fit, newdata = data.frame(X = x_test), type = "response"),
      error = function(e) rep(mean(r_train), length(x_test))
    )
    e_train <- tryCatch(
      predict(fit, newdata = data.frame(X = x_cutoff_train), type = "response"),
      error = function(e) rep(mean(r_train), length(x_cutoff_train))
    )
  }
  list(e_test = clamp_e(e_test), e_train = clamp_e(e_train))
}

make_folds <- function(n, K = 5) sample(rep(seq_len(K), length.out = n))

l2_method_grid <- function() {
  out <- expand.grid(
    power_name = names(A_POWERS),
    rho = RHO_GRID,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  out$a <- as.numeric(A_POWERS[out$power_name])
  out$label <- sprintf("L2GridTS_%s_tau%.2f_rho%.2f", out$power_name, TAU_MAX, out$rho)
  out
}

L2_METHODS <- l2_method_grid()

choose_l2grid_ec <- function(e_train, T_train, tau_max, rho_grid, a = 1, n_grid = 140) {
  e_train <- clamp_e(e_train)
  T_train <- as.numeric(T_train)
  ok <- is.finite(e_train) & is.finite(T_train)
  e <- e_train[ok]
  T <- T_train[ok]
  ntr <- length(e)

  if (ntr < 10) {
    return(data.frame(
      rho = rho_grid, ec = EPS_E, M = 0, active = 0, safe_ok = FALSE
    ))
  }

  tau_allowed <- tau_max + 1 / ntr
  probs <- unique(pmin(pmax(seq(0, tau_max, length.out = n_grid), 0), 1))
  cand <- unique(as.numeric(stats::quantile(e, probs = probs, type = 8, na.rm = TRUE)))
  cand <- sort(unique(pmax(c(EPS_E, cand), EPS_E)))

  active <- numeric(length(cand))
  M <- numeric(length(cand))
  for (jj in seq_along(cand)) {
    h <- h_stabilizer(e, cand[jj], a = a)
    active[jj] <- mean(h < 0.999, na.rm = TRUE)
    M[jj] <- sqrt(mean(((1 - h)^2) * (T^2), na.rm = TRUE))
  }

  ans <- lapply(rho_grid, function(rho) {
    feasible <- which(active <= tau_allowed & M <= rho)
    if (length(feasible) == 0) {
      data.frame(rho = rho, ec = EPS_E, M = 0, active = 0, safe_ok = FALSE)
    } else {
      best <- feasible[which.max(cand[feasible])]
      data.frame(rho = rho, ec = cand[best], M = M[best], active = active[best], safe_ok = TRUE)
    }
  })
  bind_rows(ans)
}

crossfit_nuisance_and_cutoffs <- function(X, R, Y, scenario, K = 5) {
  n <- length(X)
  folds <- make_folds(n, K)

  ehat <- rep(NA_real_, n)
  mhat <- rep(NA_real_, n)
  rhat <- rep(NA_real_, n)

  ec_quant <- rep(NA_real_, n)
  n_l2 <- nrow(L2_METHODS)
  ec_l2grid <- matrix(NA_real_, n, n_l2, dimnames = list(NULL, L2_METHODS$label))
  l2_M <- matrix(NA_real_, n, n_l2, dimnames = list(NULL, L2_METHODS$label))
  l2_active_train <- matrix(NA_real_, n, n_l2, dimnames = list(NULL, L2_METHODS$label))
  l2_ok <- matrix(NA_real_, n, n_l2, dimnames = list(NULL, L2_METHODS$label))

  tau_shrink <- tau_shrink_fun(n)

  for (k in seq_len(K)) {
    test <- which(folds == k)
    train <- which(folds != k)

    ## Propensity model fitted on training fold.
    ep <- fit_predict_e(X[train], R[train], X[test], X[train])
    e_test <- ep$e_test
    e_train <- ep$e_train
    ehat[test] <- e_test

    ## Outcome regression among observed training cases.
    obs_train <- train[R[train] == 1]
    if (length(obs_train) < 5) {
      m_const <- mean(Y[obs_train], na.rm = TRUE)
      if (!is.finite(m_const)) m_const <- mean(Y[R == 1], na.rm = TRUE)
      if (!is.finite(m_const)) m_const <- 0
      m_train <- rep(m_const, length(train))
      m_test  <- rep(m_const, length(test))
    } else {
      Bm_obs <- basis_m(X[obs_train], scenario)
      coef_m <- ridge_coef(Bm_obs, Y[obs_train])
      m_train <- as.numeric(basis_m(X[train], scenario) %*% coef_m)
      m_test  <- as.numeric(basis_m(X[test], scenario) %*% coef_m)
    }
    mhat[test] <- m_test

    ## Residual-centering regression among observed training cases.
    if (length(obs_train) < 5) {
      r_train <- rep(0, length(train))
      r_test  <- rep(0, length(test))
    } else {
      obs_pos <- match(obs_train, train)
      resid_obs <- Y[obs_train] - m_train[obs_pos]
      coef_r <- ridge_coef(basis_r(X[obs_train]), resid_obs)
      r_train <- as.numeric(basis_r(X[train]) %*% coef_r)
      r_test  <- as.numeric(basis_r(X[test]) %*% coef_r)
    }
    rhat[test] <- r_test

    ## Shrinking quantile cutoff benchmark.
    ec_q <- as.numeric(stats::quantile(e_train, probs = tau_shrink, type = 8, na.rm = TRUE))
    if (!is.finite(ec_q)) ec_q <- min(e_train, na.rm = TRUE)
    ec_q <- max(ec_q, EPS_E)
    ec_quant[test] <- ec_q

    ## L2-grid cutoffs selected on training fold.
    Z_train <- R[train] / e_train * (Y[train] - m_train)
    T_train <- Z_train - r_train

    for (pname in names(A_POWERS)) {
      a <- A_POWERS[[pname]]
      ch <- choose_l2grid_ec(e_train, T_train, tau_max = TAU_MAX, rho_grid = RHO_GRID, a = a)
      for (rr in seq_len(nrow(ch))) {
        lab <- sprintf("L2GridTS_%s_tau%.2f_rho%.2f", pname, TAU_MAX, ch$rho[rr])
        col <- match(lab, colnames(ec_l2grid))
        ec_l2grid[test, col] <- ch$ec[rr]
        l2_M[test, col] <- ch$M[rr]
        l2_active_train[test, col] <- ch$active[rr]
        l2_ok[test, col] <- as.numeric(ch$safe_ok[rr])
      }
    }
  }

  list(
    ehat = clamp_e(ehat),
    mhat = as.numeric(mhat),
    rhat = as.numeric(rhat),
    ec_quant = as.numeric(ec_quant),
    ec_l2grid = ec_l2grid,
    l2_M = l2_M,
    l2_active_train = l2_active_train,
    l2_ok = l2_ok,
    tau_shrink = tau_shrink,
    tau_max = TAU_MAX
  )
}

summarize_phi <- function(phi, est, psi0, n, h = NULL, ec = NULL,
                          est_aipw = NULL, aux = list()) {
  se <- stats::sd(phi, na.rm = TRUE) / sqrt(n)
  out <- data.frame(
    est = est,
    err = est - psi0,
    se = se,
    cover = as.numeric(is.finite(se) && abs(est - psi0) <= 1.96 * se),
    q95_abs_phi = as.numeric(stats::quantile(abs(phi), 0.95, na.rm = TRUE, type = 8)),
    q99_abs_phi = as.numeric(stats::quantile(abs(phi), 0.99, na.rm = TRUE, type = 8)),
    max_abs_phi = max(abs(phi), na.rm = TRUE),
    active_frac = if (is.null(h)) NA_real_ else mean(h < 0.999, na.rm = TRUE),
    h_l2 = if (is.null(h)) NA_real_ else sqrt(mean((h - 1)^2, na.rm = TRUE)),
    ec_mean = if (is.null(ec)) NA_real_ else mean(ec, na.rm = TRUE),
    delta_sqrt_n_from_aipw = if (is.null(est_aipw)) NA_real_ else sqrt(n) * (est - est_aipw)
  )
  if (length(aux) > 0) {
    for (nm in names(aux)) out[[nm]] <- aux[[nm]]
  }
  out
}

row_with_method <- function(method, stats_df) {
  cbind(data.frame(method = method, stringsAsFactors = FALSE), stats_df)
}

one_rep <- function(rep_id, n_samp, kappa, d_effect, scenario) {
  seed <- as.integer(
    1000000 +
      rep_id +
      10000 * match(n_samp, NSAMP) +
      100000 * match(kappa, KAPPAS) +
      1000000 * match(d_effect, D_EFFECTS) +
      10000000 * match(scenario, SCENARIOS)
  )
  set.seed(seed)

  X <- rnorm(n_samp)
  e0 <- e_fun(X, kappa)
  R <- rbinom(n_samp, 1, e0)
  mu <- mu_fun(X, d_effect)
  Y <- mu + rnorm(n_samp, sd = SIGMA)
  psi0 <- psi_true(d_effect)

  nuis <- crossfit_nuisance_and_cutoffs(X, R, Y, scenario, K = K_FOLD)
  ehat <- nuis$ehat
  mhat <- nuis$mhat
  rhat <- nuis$rhat

  ## Ordinary estimators.
  phi_or <- mhat
  est_or <- mean(phi_or)

  phi_aipw <- mhat + R / ehat * (Y - mhat)
  est_aipw <- mean(phi_aipw)

  e_trunc <- pmax(ehat, TRIM_E)
  phi_trunc <- mhat + R / e_trunc * (Y - mhat)
  est_trunc <- mean(phi_trunc)

  rows <- list()
  rows[["OR"]] <- row_with_method("OR", summarize_phi(phi_or, est_or, psi0, n_samp))
  rows[["AIPW"]] <- row_with_method("AIPW", summarize_phi(phi_aipw, est_aipw, psi0, n_samp, est_aipw = est_aipw))
  rows[["Trunc"]] <- row_with_method(
    sprintf("TruncAIPW_ec%.2f", TRIM_E),
    summarize_phi(phi_trunc, est_trunc, psi0, n_samp, est_aipw = est_aipw)
  )

  ## Fixed, shrinking-quantile, and L2-grid TS-AIPW.
  for (pname in names(A_POWERS)) {
    a <- A_POWERS[[pname]]

    h_fixed <- h_stabilizer(ehat, FIXED_EC_TS, a = a)
    phi_fixed <- mhat + h_fixed * R / ehat * (Y - mhat) + (1 - h_fixed) * rhat
    est_fixed <- mean(phi_fixed)
    lab_fixed <- sprintf("FixedTS_%s_ec%.2f", pname, FIXED_EC_TS)
    rows[[lab_fixed]] <- row_with_method(
      lab_fixed,
      summarize_phi(phi_fixed, est_fixed, psi0, n_samp,
                    h = h_fixed, ec = rep(FIXED_EC_TS, n_samp), est_aipw = est_aipw)
    )

    h_quant <- h_stabilizer(ehat, nuis$ec_quant, a = a)
    phi_quant <- mhat + h_quant * R / ehat * (Y - mhat) + (1 - h_quant) * rhat
    est_quant <- mean(phi_quant)
    lab_quant <- sprintf("QuantileTS_%s_shrink", pname)
    rows[[lab_quant]] <- row_with_method(
      lab_quant,
      summarize_phi(phi_quant, est_quant, psi0, n_samp,
                    h = h_quant, ec = nuis$ec_quant, est_aipw = est_aipw)
    )
  }

  for (jj in seq_len(nrow(L2_METHODS))) {
    lab <- L2_METHODS$label[jj]
    a <- L2_METHODS$a[jj]
    ec_grid <- nuis$ec_l2grid[, lab]
    h_grid <- h_stabilizer(ehat, ec_grid, a = a)
    phi_grid <- mhat + h_grid * R / ehat * (Y - mhat) + (1 - h_grid) * rhat
    est_grid <- mean(phi_grid)
    rows[[lab]] <- row_with_method(
      lab,
      summarize_phi(
        phi_grid, est_grid, psi0, n_samp,
        h = h_grid, ec = ec_grid, est_aipw = est_aipw,
        aux = list(
          l2_tau_max = TAU_MAX,
          l2_rho = L2_METHODS$rho[jj],
          l2_M_train = mean(nuis$l2_M[, lab], na.rm = TRUE),
          l2_active_train = mean(nuis$l2_active_train[, lab], na.rm = TRUE),
          l2_ok_rate = mean(nuis$l2_ok[, lab], na.rm = TRUE)
        )
      )
    )
  }

  out <- dplyr::bind_rows(rows) %>%
    dplyr::mutate(
      rep_id = rep_id,
      n_samp = n_samp,
      kappa = kappa,
      d_effect = d_effect,
      scenario = scenario,
      psi0 = psi0,
      tau_shrink = nuis$tau_shrink,
      tau_max = nuis$tau_max,
      response_rate = mean(R),
      e_min = min(ehat),
      e_q01 = as.numeric(stats::quantile(ehat, 0.01, na.rm = TRUE, type = 8)),
      e_q05 = as.numeric(stats::quantile(ehat, 0.05, na.rm = TRUE, type = 8))
    ) %>%
    dplyr::relocate(rep_id, scenario, n_samp, kappa, d_effect, method)

  out
}

## ==========================================================
## 3. Run simulations in parallel with checkpointing
## ==========================================================

jobs <- expand.grid(
  rep_id = seq_len(N_REP),
  scenario = SCENARIOS,
  n_samp = NSAMP,
  kappa = KAPPAS,
  d_effect = D_EFFECTS,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

jobs <- dplyr::as_tibble(jobs) %>%
  dplyr::arrange(scenario, kappa, d_effect, n_samp, rep_id)

cat("Number of jobs:", nrow(jobs), "\n")
cat("Number of method rows expected:", nrow(jobs) * (4 + 2 * length(A_POWERS) + nrow(L2_METHODS)), "\n\n")

run_job <- function(j) {
  tryCatch(
    one_rep(
      rep_id = jobs$rep_id[j],
      n_samp = jobs$n_samp[j],
      kappa = jobs$kappa[j],
      d_effect = jobs$d_effect[j],
      scenario = jobs$scenario[j]
    ),
    error = function(e) {
      data.frame(
        rep_id = jobs$rep_id[j],
        scenario = jobs$scenario[j],
        n_samp = jobs$n_samp[j],
        kappa = jobs$kappa[j],
        d_effect = jobs$d_effect[j],
        method = "ERROR",
        error_message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    }
  )
}

CHUNK_SIZE <- if (QUICK_TEST) 50 else 300
chunks <- split(seq_len(nrow(jobs)), ceiling(seq_along(seq_len(nrow(jobs))) / CHUNK_SIZE))

detail_all <- list()
t0 <- Sys.time()

for (cc in seq_along(chunks)) {
  idx <- chunks[[cc]]
  cat("\n============================================================\n")
  cat("Chunk", cc, "of", length(chunks), "jobs", min(idx), "--", max(idx), "\n")
  cat("Time:", format(Sys.time()), "\n")
  cat("============================================================\n")

  if (.Platform$OS.type == "unix" && N_CORES > 1) {
    res <- parallel::mclapply(idx, run_job, mc.cores = N_CORES, mc.preschedule = FALSE)
  } else {
    cl <- parallel::makeCluster(N_CORES)
    parallel::clusterExport(cl, varlist = setdiff(ls(envir = .GlobalEnv), c("cl")), envir = .GlobalEnv)
    res <- parallel::parLapply(cl, idx, run_job)
    parallel::stopCluster(cl)
  }

  detail_all[[cc]] <- dplyr::bind_rows(res)
  detail_tmp <- dplyr::bind_rows(detail_all)
  readr::write_csv(detail_tmp, file.path(OUTDIR, "adaptive_TS_L2grid_detail_checkpoint.csv"))

  cat("Rows saved so far:", nrow(detail_tmp), "\n")
  cat("Elapsed minutes:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
}

detail <- dplyr::bind_rows(detail_all)
readr::write_csv(detail, file.path(OUTDIR, "adaptive_TS_L2grid_detail.csv"))

if (any(detail$method == "ERROR", na.rm = TRUE)) {
  cat("\nErrors detected. See rows with method == ERROR in adaptive_TS_L2grid_detail.csv\n")
  print(detail %>% dplyr::filter(method == "ERROR") %>% head(10))
}

detail_ok <- detail %>%
  dplyr::filter(method != "ERROR") %>%
  dplyr::mutate(method = as.character(method), scenario = as.character(scenario))

summary <- detail_ok %>%
  dplyr::group_by(scenario, n_samp, kappa, d_effect, method) %>%
  dplyr::summarise(
    n_rep = dplyr::n(),
    est_mean = mean(est, na.rm = TRUE),
    mean_bias = mean(est - psi0, na.rm = TRUE),
    emp_sd = sd(est, na.rm = TRUE),
    rmse = sqrt(mean((est - psi0)^2, na.rm = TRUE)),
    mean_se = mean(se, na.rm = TRUE),
    coverage = mean(cover, na.rm = TRUE),
    q95_abs_phi_mean = mean(q95_abs_phi, na.rm = TRUE),
    q99_abs_phi_mean = mean(q99_abs_phi, na.rm = TRUE),
    max_abs_phi_median = median(max_abs_phi, na.rm = TRUE),
    active_frac_mean = mean(active_frac, na.rm = TRUE),
    h_l2_mean = mean(h_l2, na.rm = TRUE),
    ec_mean = mean(ec_mean, na.rm = TRUE),
    tau_shrink = mean(tau_shrink, na.rm = TRUE),
    tau_max = mean(tau_max, na.rm = TRUE),
    l2_rho = mean(l2_rho, na.rm = TRUE),
    delta_mean = mean(delta_sqrt_n_from_aipw, na.rm = TRUE),
    delta_sd = sd(delta_sqrt_n_from_aipw, na.rm = TRUE),
    delta_mad = median(abs(delta_sqrt_n_from_aipw), na.rm = TRUE),
    delta_q90 = as.numeric(stats::quantile(abs(delta_sqrt_n_from_aipw), 0.90, na.rm = TRUE, type = 8)),
    delta_q95 = as.numeric(stats::quantile(abs(delta_sqrt_n_from_aipw), 0.95, na.rm = TRUE, type = 8)),
    delta_q99 = as.numeric(stats::quantile(abs(delta_sqrt_n_from_aipw), 0.99, na.rm = TRUE, type = 8)),
    l2_M_train = mean(l2_M_train, na.rm = TRUE),
    l2_active_train = mean(l2_active_train, na.rm = TRUE),
    l2_ok_rate = mean(l2_ok_rate, na.rm = TRUE),
    response_rate = mean(response_rate, na.rm = TRUE),
    e_q01 = mean(e_q01, na.rm = TRUE),
    e_q05 = mean(e_q05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(scenario, kappa, d_effect, n_samp, method)

readr::write_csv(summary, file.path(OUTDIR, "adaptive_TS_L2grid_summary.csv"))

cat("\n============================================================\n")
cat("Finished adaptive / L2-grid TS-AIPW simulation\n")
cat("Elapsed total minutes:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("Detail :", file.path(OUTDIR, "adaptive_TS_L2grid_detail.csv"), "\n")
cat("Summary:", file.path(OUTDIR, "adaptive_TS_L2grid_summary.csv"), "\n")
cat("============================================================\n\n")

## ==========================================================
## 4. Manuscript-ready reduced tables
## ==========================================================

main_methods <- c(
  "AIPW",
  sprintf("TruncAIPW_ec%.2f", TRIM_E),
  sprintf("FixedTS_a1_ec%.2f", FIXED_EC_TS),
  "QuantileTS_a1_shrink",
  L2_METHODS$label
)

method_label <- function(x) {
  out <- x
  out[x == "AIPW"] <- "AIPW"
  out[x == sprintf("TruncAIPW_ec%.2f", TRIM_E)] <- "Truncated AIPW"
  out[x == sprintf("FixedTS_a1_ec%.2f", FIXED_EC_TS)] <- "Fixed TS-AIPW"
  out[x == "QuantileTS_a1_shrink"] <- "Shrinking quantile TS"
  for (jj in seq_len(nrow(L2_METHODS))) {
    lab <- L2_METHODS$label[jj]
    out[x == lab] <- sprintf("L2-grid TS (rho=%.2f)", L2_METHODS$rho[jj])
  }
  out
}

method_order_tbl <- tibble(method = main_methods, method_order = seq_along(main_methods))

table_finite <- summary %>%
  dplyr::filter(scenario == "coarse_m", method %in% main_methods) %>%
  dplyr::left_join(method_order_tbl, by = "method") %>%
  dplyr::mutate(method_label = method_label(method)) %>%
  dplyr::arrange(kappa, d_effect, n_samp, method_order) %>%
  dplyr::select(
    n_samp, kappa, d_effect, method = method_label,
    mean_bias, emp_sd, rmse, coverage,
    q99_abs_phi_mean, max_abs_phi_median,
    active_frac_mean, h_l2_mean,
    delta_sd, delta_mad, delta_q95, delta_q99,
    tau_max, l2_rho, l2_M_train, l2_active_train, l2_ok_rate
  )

readr::write_csv(table_finite, file.path(OUTDIR, "table_finite_L2grid.csv"))

table_eff <- summary %>%
  dplyr::filter(
    scenario == "correct_m",
    method %in% c("AIPW", "QuantileTS_a1_shrink", sprintf("FixedTS_a1_ec%.2f", FIXED_EC_TS), L2_METHODS$label)
  ) %>%
  dplyr::left_join(method_order_tbl, by = "method") %>%
  dplyr::mutate(method_label = method_label(method)) %>%
  dplyr::arrange(d_effect, n_samp, method_order) %>%
  dplyr::select(
    n_samp, kappa, d_effect, method = method_label,
    mean_bias, emp_sd, mean_se, rmse, coverage,
    active_frac_mean, h_l2_mean,
    delta_sd, delta_mad, delta_q95, delta_q99,
    tau_max, l2_rho, l2_M_train, l2_active_train, l2_ok_rate
  )

readr::write_csv(table_eff, file.path(OUTDIR, "table_efficiency_L2grid.csv"))

cat("\n===== Finite-sample L2-grid table: coarse_m =====\n")
print(tibble::as_tibble(table_finite), n = 120, width = Inf)

cat("\n===== Efficiency / recovery table: correct_m =====\n")
print(tibble::as_tibble(table_eff), n = 120, width = Inf)

## ==========================================================
## 5. Figures, with empty-data guards
## ==========================================================

theme_paper <- function(base_size = 13) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey70")
    )
}

plot_df <- summary %>%
  dplyr::filter(method %in% main_methods) %>%
  dplyr::mutate(method_label = method_label(method))

save_plot_if_data <- function(df, filename, plot_fun, width = 8, height = 4.8) {
  if (nrow(df) == 0) {
    cat("Skipping", filename, "because there are no rows.\n")
    return(invisible(NULL))
  }
  p <- plot_fun(df)
  ggplot2::ggsave(file.path(OUTDIR, filename), p, width = width, height = height, dpi = 300)
  cat("Saved", file.path(OUTDIR, filename), "\n")
}

## Figure 1: RMSE under approximate residual recovery.
df_rmse <- plot_df %>%
  dplyr::filter(
    scenario == "coarse_m",
    method_label %in% c("AIPW", "Truncated AIPW", "Fixed TS-AIPW",
                        "Shrinking quantile TS", "L2-grid TS (rho=0.05)",
                        "L2-grid TS (rho=0.10)", "L2-grid TS (rho=0.20)",
                        "L2-grid TS (rho=0.40)")
  )
save_plot_if_data(df_rmse, "fig_rmse_coarse_L2grid.png", function(df) {
  ggplot2::ggplot(df, ggplot2::aes(x = n_samp, y = rmse,
                                   group = method_label, linetype = method_label,
                                   shape = method_label)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_log10(breaks = NSAMP) +
    ggplot2::labs(
      x = "Sample size n",
      y = "RMSE",
      linetype = "Method",
      shape = "Method",
      title = "Finite-sample stabilization with L2-grid adaptive cutoff"
    ) +
    ggplot2::facet_wrap(~ d_effect, scales = "free_y",
                        labeller = ggplot2::label_bquote(d == .(d_effect))) +
    theme_paper()
}, width = 10, height = 5)

## Figure 2: tail quantile.
df_tail <- df_rmse
save_plot_if_data(df_tail, "fig_tail_q99_coarse_L2grid.png", function(df) {
  ggplot2::ggplot(df, ggplot2::aes(x = n_samp, y = q99_abs_phi_mean,
                                   group = method_label, linetype = method_label,
                                   shape = method_label)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_log10(breaks = NSAMP) +
    ggplot2::labs(
      x = "Sample size n",
      y = "mean q.99(|pseudo-outcome|)",
      linetype = "Method",
      shape = "Method",
      title = "Pseudo-outcome tail diagnostic"
    ) +
    ggplot2::facet_wrap(~ d_effect, scales = "free_y",
                        labeller = ggplot2::label_bquote(d == .(d_effect))) +
    theme_paper()
}, width = 10, height = 5)

## Figure 3: delta diagnostic against AIPW.
df_delta <- plot_df %>%
  dplyr::filter(
    method_label %in% c("Fixed TS-AIPW", "Shrinking quantile TS",
                        "L2-grid TS (rho=0.05)", "L2-grid TS (rho=0.10)",
                        "L2-grid TS (rho=0.20)", "L2-grid TS (rho=0.40)"),
    d_effect %in% max(D_EFFECTS)
  )
save_plot_if_data(df_delta, "fig_delta_q95_L2grid.png", function(df) {
  ggplot2::ggplot(df, ggplot2::aes(x = n_samp, y = delta_q95,
                                   group = method_label, linetype = method_label,
                                   shape = method_label)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_log10(breaks = NSAMP) +
    ggplot2::labs(
      x = "Sample size n",
      y = "q.95 |sqrt(n)(psi_TS - psi_AIPW)|",
      linetype = "Method",
      shape = "Method",
      title = "Deviation from ordinary AIPW"
    ) +
    ggplot2::facet_wrap(~ scenario, scales = "free_y") +
    theme_paper()
}, width = 10, height = 5)

## Figure 4: selected training deleted-residual L2 mass.
df_l2 <- plot_df %>%
  dplyr::filter(stringr::str_detect(method_label, "L2-grid TS"))
save_plot_if_data(df_l2, "fig_L2grid_train_mass.png", function(df) {
  ggplot2::ggplot(df, ggplot2::aes(x = n_samp, y = l2_M_train,
                                   group = method_label, linetype = method_label,
                                   shape = method_label)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_log10(breaks = NSAMP) +
    ggplot2::labs(
      x = "Sample size n",
      y = "training deleted-residual L2 mass",
      linetype = "Method",
      shape = "Method",
      title = "L2-grid cutoff diagnostic"
    ) +
    ggplot2::facet_grid(scenario ~ d_effect, scales = "free_y") +
    theme_paper()
}, width = 11, height = 6)

cat("\nAll done.\n")
