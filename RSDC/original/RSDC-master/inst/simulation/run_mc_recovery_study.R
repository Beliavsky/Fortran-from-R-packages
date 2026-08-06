# =============================================================================
# Monte Carlo Study - RSDC Package (1.7-0, global task queue)
# =============================================================================
#
# Five cases:
#   Case 1 - const, N=1, K=2
#   Case 2 - noX,   N=2, K=2
#   Case 3 - tvtp,  N=2, K=2
#   Case 4 - noX,   N=3, K in {2, 3, 5}   (K = 5: the 1.7-0 search at work)
#   Case 5 - tvtp,  N=3, K=2              (well-identified design)
#
# Design: M = 100 replications; T in {500, 1000, 2000} (cases 1-4),
# {1000, 2000} (case 5).
#
# Relation to the previous study (pre-1.7-0 optimizer):
#   * The DGPs, true parameters and BASE SEEDS are unchanged, so every cell
#     that existed before re-estimates the *identical simulated data sets* -
#     differences in bias/RMSE/usable rates isolate the effect of the
#     reparameterized global search (partial correlations + bounded softmax
#     + top-3 refinement) introduced in 1.7-0.
#   * rsdc_simulate() and the estimator share the same transition links
#     (logistic N=2 / softmax with N-th reference logit for N>=3), the same
#     X[t] timing and the same t=1 prior, so simulator and likelihood tell
#     the same story; nothing in the DGPs needed to change.
#   * compute_se = FALSE (SEs play no role in bias/RMSE; estimates identical).
#   * Two estimation arms: "default" (single search, out-of-the-box, as the
#     previous study) for all cells, plus "ns4" (control$n_starts = 4, the
#     documented protocol for multimodal surfaces) for the N = 3 cases 4-5.
#     Same data in both arms: the delta isolates the multi-start's value.
#
# Parallelization: ALL (case, K, T, m) replications are flattened into one
# global task queue mapped over the workers by a single parallel::mclapply
# with dynamic scheduling (mc.preschedule = FALSE) - long tvtp N=3 tasks
# interleave with fast const ones, no per-cell barriers, no idle tails.
# Each task self-seeds (set.seed(base_seed + m)): results are bit-identical
# for any worker count. Backend: forked mclapply on macOS/Linux, a PSOCK
# cluster on Windows.
#
# Outputs: PNG plots -> inst/simulation/plots/
#          CSV summaries -> inst/simulation/results/
#          RDS bundle -> inst/simulation/results/monte_carlo_results.rds
#
# Env overrides: MC_CORES (default: available cores - 2), MC_M (default 100),
#                MC_SMOKE=1 (tiny grids, syntax/pipeline check only).
#
# Source from the package root: source("inst/simulation/run_mc_recovery_study.R")
# =============================================================================

# How the package is loaded here is replayed verbatim on PSOCK workers below
# (forked workers inherit it, PSOCK workers do not), so a development
# load_all() session parallelises the same way an installed package does.
LOAD_RSDC <- if (requireNamespace("RSDC", quietly = TRUE)) {
  function() library("RSDC")
} else if (requireNamespace("devtools", quietly = TRUE)) {
  root <- tryCatch(rprojroot::find_root(rprojroot::has_file("DESCRIPTION")),
                   error = function(e) getwd())
  function() devtools::load_all(root)
} else {
  stop("Install RSDC (or 'devtools' for development) before running this script.")
}
LOAD_RSDC()
stopifnot(utils::packageVersion("RSDC") >= "1.7.0")

# -- Global configuration ------------------------------------------------------
SMOKE    <- identical(Sys.getenv("MC_SMOKE"), "1")
M        <- {
  env <- suppressWarnings(as.integer(Sys.getenv("MC_M", "")))
  if (!is.na(env) && env >= 1L) env else if (SMOKE) 2L else 100L
}
T_grid   <- if (SMOKE) 300L else c(500L, 1000L, 2000L)   # cases 1-4
T5_grid  <- if (SMOKE) 300L else c(1000L, 2000L)         # case 5 (local; never
                                                         # clobbers T_grid)
K4_grid  <- if (SMOKE) 2L   else c(2L, 3L, 5L)           # case-4 K sweep
MC_CORES <- {
  env <- suppressWarnings(as.integer(Sys.getenv("MC_CORES", "")))
  if (!is.na(env) && env >= 1L) env
  else {
    # detectCores() may return NA on some platforms; fall back to serial.
    dc <- parallel::detectCores()
    if (is.na(dc) || !is.finite(dc)) 1L else max(1L, dc - 2L)
  }
}
CTRL <- list(compute_se = FALSE)     # SEs play no role in bias/RMSE
# Second arm for the multimodal N = 3 cases (4 and 5): the documented
# protocol for multimodal surfaces (n_starts = 4 independent searches, best
# kept). Same seeds -> same simulated data; only the estimation protocol
# differs, so "default vs protocol" isolates the value of the multi-start.
CTRL_NS4 <- list(compute_se = FALSE, n_starts = 4L)

# Pin BLAS to one thread under the M-way fork (no-op if RhpcBLASctl absent).
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1L)

plots_dir   <- file.path("inst", "simulation", "plots")
results_dir <- file.path("inst", "simulation", "results")
dir.create(plots_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
saved_files <- character(0)

# -- True parameters (identical to the previous study) -------------------------
C1_RHO <- 0.5; C1_K <- 2L

C2_P11 <- 0.9; C2_P22 <- 0.8; C2_RHO1 <- 0.2; C2_RHO2 <- 0.7
C2_N <- 2L; C2_K <- 2L
C2_BETA <- matrix(c(qlogis(C2_P11), qlogis(C2_P22)), nrow = 2L, ncol = 1L)

C3_BETA <- matrix(c(qlogis(0.70), -0.8,
                    qlogis(0.70),  0.8), nrow = 2L, byrow = TRUE)
C3_RHO1 <- 0.2; C3_RHO2 <- 0.7
C3_N <- 2L; C3_K <- 2L; C3_PHI <- 0.7
C3_P11 <- plogis(C3_BETA[1L, 1L]); C3_P22 <- plogis(C3_BETA[2L, 1L])

# Case 4 (noX, N=3): ascending-rho regimes, equicorrelation per regime, the
# intercept-only softmax beta reproduces C4_P exactly (reference = regime 3).
C4_N <- 3L
C4_P <- matrix(c(0.90, 0.07, 0.03,
                 0.08, 0.85, 0.07,
                 0.05, 0.15, 0.80), nrow = 3L, byrow = TRUE)
C4_RHO  <- c(0.2, 0.5, 0.8)
C4_BETA <- t(apply(C4_P, 1L, function(pr) c(log(pr[1] / pr[3]), log(pr[2] / pr[3]))))

# Case 5 (tvtp, N=3, K=2): well-identified design - well-separated
# correlations, persistent regimes, AR(1) covariate.
C5_N <- 3L; C5_K <- 2L; C5_PHI <- 0.7
C5_RHO <- c(-0.6, 0.1, 0.8)
C5_P0  <- matrix(c(0.90, 0.07, 0.03,
                   0.05, 0.90, 0.05,
                   0.03, 0.07, 0.90), nrow = 3L, byrow = TRUE)
C5_SLOPE <- matrix(c( 0.5, 0.0,
                      0.0, 0.5,
                     -0.5, 0.0), nrow = 3L, byrow = TRUE)
# beta = N x (N-1)*p (p = 2), per row (int1, slope1, int2, slope2), ref = regime 3
C5_BETA <- t(vapply(1:3, function(i) {
  ints <- c(log(C5_P0[i, 1] / C5_P0[i, 3]), log(C5_P0[i, 2] / C5_P0[i, 3]))
  c(ints[1], C5_SLOPE[i, 1], ints[2], C5_SLOPE[i, 2])
}, numeric(4L)))

# -- Shared helpers ------------------------------------------------------------
make_corr <- function(K, rho_vals) {
  R <- diag(K)
  R[lower.tri(R)] <- rho_vals
  R[upper.tri(R)] <- t(R)[upper.tri(R)]
  R
}

make_summary <- function(est_mat, true_vals) {
  means <- colMeans(est_mat)
  data.frame(param = names(true_vals), true = unname(true_vals),
             mean_hat = means, bias = means - true_vals,
             rmse = sqrt(colMeans(sweep(est_mat, 2L, true_vals)^2)),
             sd = apply(est_mat, 2L, sd), row.names = NULL)
}

print_summary <- function(smry) { cat("\n"); print(smry, digits = 4L, row.names = FALSE); invisible(smry) }

save_csv <- function(df, fname) {
  path <- file.path(results_dir, fname)
  utils::write.csv(df, path, row.names = FALSE)
  saved_files <<- c(saved_files, path)
  invisible(path)
}

# Save a PNG; also draw on the active device when running interactively.
draw_and_save <- function(plot_fun, fname, width = 1400, height = 900, res = 130) {
  if (interactive()) plot_fun()
  path <- file.path(plots_dir, fname)
  grDevices::png(path, width = width, height = height, res = res)
  plot_fun()
  grDevices::dev.off()
  saved_files <<- c(saved_files, path)
  invisible(path)
}

boxplot_grid <- function(est_mat, true_vals, cell_lab, mfrow, col) {
  par(mfrow = mfrow)
  for (nm in colnames(est_mat)) {
    boxplot(est_mat[, nm], main = sprintf("%s | %s", nm, cell_lab),
            col = col, ylab = nm)
    abline(h = true_vals[[nm]], col = "red", lwd = 2)
  }
  par(mfrow = c(1, 1))
}

# AR(1) covariate (caller controls the RNG)
ar1_X <- function(T, phi) {
  x <- numeric(T); x[1] <- rnorm(1L, sd = 1 / sqrt(1 - phi^2))
  for (t in 2:T) x[t] <- phi * x[t - 1L] + rnorm(1L)
  cbind(1, as.numeric(scale(x)))
}

# -- Per-replication workers (one per case) ------------------------------------
# Each returns a one-row data.frame; crashes are caught and replaced by the
# case's FAIL row during aggregation. Base-seed formulas are UNCHANGED from
# the previous study, so identical data sets are re-estimated.

FAIL <- list(
  case1 = data.frame(converged = FALSE, rho_hat = NA_real_),
  case2 = data.frame(converged = FALSE, p11_hat = NA_real_, p22_hat = NA_real_,
                     rho1_hat = NA_real_, rho2_hat = NA_real_),
  case3 = data.frame(converged = FALSE, label_ok = FALSE, bound_hit = NA,
                     p11_hat = NA_real_, p22_hat = NA_real_,
                     b11_hat = NA_real_, b12_hat = NA_real_,
                     b21_hat = NA_real_, b22_hat = NA_real_,
                     rho1_hat = NA_real_, rho2_hat = NA_real_),
  case4 = data.frame(converged = FALSE,
                     p11_hat = NA_real_, p22_hat = NA_real_, p33_hat = NA_real_,
                     rho1_hat = NA_real_, rho2_hat = NA_real_, rho3_hat = NA_real_),
  case5 = data.frame(converged = FALSE, label_ok = FALSE, bound_hit = NA,
                     p11_hat = NA_real_, p22_hat = NA_real_, p33_hat = NA_real_,
                     b11_hat = NA_real_, b12_hat = NA_real_, b13_hat = NA_real_,
                     b14_hat = NA_real_,
                     b21_hat = NA_real_, b22_hat = NA_real_, b23_hat = NA_real_,
                     b24_hat = NA_real_,
                     b31_hat = NA_real_, b32_hat = NA_real_, b33_hat = NA_real_,
                     b34_hat = NA_real_,
                     rho1_hat = NA_real_, rho2_hat = NA_real_, rho3_hat = NA_real_)
)

run_case1 <- function(m, TT, KK, ctrl = CTRL) {
  set.seed(20000L + TT * 1000L + m)
  Sigma <- array(make_corr(C1_K, C1_RHO), dim = c(C1_K, C1_K, 1L))
  sim <- rsdc_simulate(n = TT, X = matrix(1, TT, 1L),
                       beta = matrix(nrow = 1L, ncol = 0L),
                       mu = matrix(0, 1L, C1_K), sigma = Sigma, N = 1L)
  fit <- rsdc_estimate("const", residuals = sim$observations, control = ctrl)
  data.frame(converged = TRUE, rho_hat = fit$correlations[1L, 1L])
}

run_case2 <- function(m, TT, KK, ctrl = CTRL) {
  set.seed(10000L + TT * 1000L + m)
  Sigma <- array(0, dim = c(C2_K, C2_K, C2_N))
  Sigma[,, 1L] <- make_corr(C2_K, C2_RHO1); Sigma[,, 2L] <- make_corr(C2_K, C2_RHO2)
  sim <- rsdc_simulate(n = TT, X = matrix(1, TT, 1L), beta = C2_BETA,
                       mu = matrix(0, C2_N, C2_K), sigma = Sigma, N = C2_N)
  fit <- rsdc_estimate("noX", residuals = sim$observations, N = C2_N, control = ctrl)
  data.frame(converged = TRUE,
             p11_hat = fit$transition_matrix[1L, 1L],
             p22_hat = fit$transition_matrix[2L, 2L],
             rho1_hat = fit$correlations[1L, 1L],
             rho2_hat = fit$correlations[2L, 1L])
}

run_case3 <- function(m, TT, KK, ctrl = CTRL) {
  set.seed(31000L + TT * 1000L + m)
  X <- ar1_X(TT, C3_PHI)
  Sigma <- array(0, dim = c(C3_K, C3_K, C3_N))
  Sigma[,, 1L] <- make_corr(C3_K, C3_RHO1); Sigma[,, 2L] <- make_corr(C3_K, C3_RHO2)
  sim <- rsdc_simulate(n = TT, X = X, beta = C3_BETA,
                       mu = matrix(0, C3_N, C3_K), sigma = Sigma, N = C3_N)
  fit <- rsdc_estimate("tvtp", residuals = sim$observations, N = C3_N, X = X,
                       control = ctrl)
  r1 <- fit$correlations[1L, 1L]; r2 <- fit$correlations[2L, 1L]
  data.frame(converged = TRUE, label_ok = r2 > r1,
             bound_hit = any(abs(fit$beta) >= 9.5),
             p11_hat = fit$transition_matrix[1L, 1L],
             p22_hat = fit$transition_matrix[2L, 2L],
             b11_hat = fit$beta[1L, 1L], b12_hat = fit$beta[1L, 2L],
             b21_hat = fit$beta[2L, 1L], b22_hat = fit$beta[2L, 2L],
             rho1_hat = r1, rho2_hat = r2)
}

run_case4 <- function(m, TT, KK, ctrl = CTRL) {
  set.seed(40000000L + KK * 3000000L + TT * 1000L + m)
  Sigma <- array(0, dim = c(KK, KK, C4_N))
  for (i in seq_len(C4_N))
    Sigma[,, i] <- make_corr(KK, rep(C4_RHO[i], KK * (KK - 1L) / 2L))
  sim <- rsdc_simulate(n = TT, X = matrix(1, TT, 1L), beta = C4_BETA,
                       mu = matrix(0, C4_N, KK), sigma = Sigma, N = C4_N)
  fit <- rsdc_estimate("noX", residuals = sim$observations, N = C4_N, control = ctrl)
  data.frame(converged = TRUE,
             p11_hat = fit$transition_matrix[1L, 1L],
             p22_hat = fit$transition_matrix[2L, 2L],
             p33_hat = fit$transition_matrix[3L, 3L],
             rho1_hat = mean(fit$correlations[1L, ]),
             rho2_hat = mean(fit$correlations[2L, ]),
             rho3_hat = mean(fit$correlations[3L, ]))
}

run_case5 <- function(m, TT, KK, ctrl = CTRL) {
  set.seed(70000000L + KK * 3000000L + TT * 1000L + m)
  X <- ar1_X(TT, C5_PHI)
  Sigma <- array(0, dim = c(KK, KK, C5_N))
  for (i in seq_len(C5_N))
    Sigma[,, i] <- make_corr(KK, rep(C5_RHO[i], KK * (KK - 1L) / 2L))
  sim <- rsdc_simulate(n = TT, X = X, beta = C5_BETA,
                       mu = matrix(0, C5_N, KK), sigma = Sigma, N = C5_N)
  fit <- rsdc_estimate("tvtp", residuals = sim$observations, N = C5_N, X = X,
                       control = ctrl)
  r <- rowMeans(fit$correlations)
  data.frame(converged = TRUE, label_ok = (r[1] < r[2]) & (r[2] < r[3]),
             bound_hit = any(abs(fit$beta) >= 9.5),
             p11_hat = fit$transition_matrix[1L, 1L],
             p22_hat = fit$transition_matrix[2L, 2L],
             p33_hat = fit$transition_matrix[3L, 3L],
             b11_hat = fit$beta[1L, 1L], b12_hat = fit$beta[1L, 2L],
             b13_hat = fit$beta[1L, 3L], b14_hat = fit$beta[1L, 4L],
             b21_hat = fit$beta[2L, 1L], b22_hat = fit$beta[2L, 2L],
             b23_hat = fit$beta[2L, 3L], b24_hat = fit$beta[2L, 4L],
             b31_hat = fit$beta[3L, 1L], b32_hat = fit$beta[3L, 2L],
             b33_hat = fit$beta[3L, 3L], b34_hat = fit$beta[3L, 4L],
             rho1_hat = r[1], rho2_hat = r[2], rho3_hat = r[3])
}

WORKERS <- list(case1 = run_case1, case2 = run_case2, case3 = run_case3,
                case4 = run_case4, case5 = run_case5)

# -- Global task queue ---------------------------------------------------------
# One row per replication across ALL cells; a single dynamic mclapply maps the
# queue over the workers. Long tasks are enqueued FIRST (largest K*T of the
# N = 3 cases) so they start early and cannot serialize the tail.
cells <- rbind(
  expand.grid(case = "case1", K = C1_K, T = T_grid,  stringsAsFactors = FALSE),
  expand.grid(case = "case2", K = C2_K, T = T_grid,  stringsAsFactors = FALSE),
  expand.grid(case = "case3", K = C3_K, T = T_grid,  stringsAsFactors = FALSE),
  expand.grid(case = "case4", K = K4_grid, T = T_grid, stringsAsFactors = FALSE),
  expand.grid(case = "case5", K = C5_K, T = T5_grid, stringsAsFactors = FALSE))
cells$arm <- "default"
# Protocol arm (n_starts = 4) for the multimodal N = 3 cases only.
ns4 <- cells[cells$case %in% c("case4", "case5"), ]
ns4$arm <- "ns4"
cells <- rbind(cells, ns4)
tasks <- merge(cells, data.frame(m = seq_len(M)))
weight <- with(tasks, ifelse(case %in% c("case4", "case5"), 10, 1) * K * T *
                      ifelse(arm == "ns4", 4, 1))
tasks  <- tasks[order(-weight), ]

cat(sprintf("Monte Carlo: RSDC %s | %d cells x M=%d -> %d tasks | cores=%d%s\n",
            as.character(utils::packageVersion("RSDC")), nrow(cells), M,
            nrow(tasks), MC_CORES, if (SMOKE) " [SMOKE]" else ""))

t0_all <- proc.time()["elapsed"]
run_task <- function(i) {
  tk <- tasks[i, ]
  out <- tryCatch(WORKERS[[tk$case]](tk$m, tk$T, tk$K,
                                     ctrl = if (tk$arm == "ns4") CTRL_NS4 else CTRL),
                  error = function(e) FAIL[[tk$case]])
  cbind(tk, out, row.names = NULL)
}
raw <- if (.Platform$OS.type == "unix") {
  parallel::mclapply(seq_len(nrow(tasks)), run_task,
                     mc.cores = MC_CORES, mc.preschedule = FALSE)
} else {
  # Windows: PSOCK cluster with load balancing. parLapplyLB returns results in
  # input order and every task self-seeds, so the output matches the forked
  # path exactly. Globals are exported first, then each worker loads RSDC the
  # same way the master did (installed package or development load_all).
  cl <- parallel::makeCluster(MC_CORES)
  tryCatch({
    parallel::clusterExport(cl, ls(envir = .GlobalEnv), envir = .GlobalEnv)
    parallel::clusterEvalQ(cl, {
      suppressMessages(LOAD_RSDC())
      # Forked workers inherit these from the master; PSOCK workers must be
      # told, or BLAS threads oversubscribe and the RNG kind may differ.
      if (requireNamespace("RhpcBLASctl", quietly = TRUE))
        RhpcBLASctl::blas_set_num_threads(1L)
      RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    })
    # chunk.size = 1L: without it parLapplyLB pre-splits into 2 * nworkers
    # static chunks, which would undo the descending-weight ordering above.
    parallel::parLapplyLB(cl, seq_len(nrow(tasks)), run_task, chunk.size = 1L)
  }, finally = parallel::stopCluster(cl))
}
elapsed_all <- proc.time()["elapsed"] - t0_all
ok_runs <- vapply(raw, is.data.frame, logical(1L))
if (any(!ok_runs))       # crashed workers (fork kill etc.): mark as failed
  raw[!ok_runs] <- lapply(which(!ok_runs), function(i)
    cbind(tasks[i, ], FAIL[[tasks$case[i]]], row.names = NULL))
cat(sprintf("All tasks done in %.1f min.\n", elapsed_all / 60))

# -- Aggregation per cell ------------------------------------------------------
TRUES <- list(
  case1 = c(rho = C1_RHO),
  case2 = c(p11 = C2_P11, p22 = C2_P22, rho1 = C2_RHO1, rho2 = C2_RHO2),
  case3 = c(p11 = C3_P11, p22 = C3_P22,
            b11 = C3_BETA[1, 1], b12 = C3_BETA[1, 2],
            b21 = C3_BETA[2, 1], b22 = C3_BETA[2, 2],
            rho1 = C3_RHO1, rho2 = C3_RHO2),
  case4 = c(p11 = C4_P[1, 1], p22 = C4_P[2, 2], p33 = C4_P[3, 3],
            rho1 = C4_RHO[1], rho2 = C4_RHO[2], rho3 = C4_RHO[3]),
  case5 = c(p11 = C5_P0[1, 1], p22 = C5_P0[2, 2], p33 = C5_P0[3, 3],
            b11 = C5_BETA[1, 1], b12 = C5_BETA[1, 2], b13 = C5_BETA[1, 3],
            b14 = C5_BETA[1, 4],
            b21 = C5_BETA[2, 1], b22 = C5_BETA[2, 2], b23 = C5_BETA[2, 3],
            b24 = C5_BETA[2, 4],
            b31 = C5_BETA[3, 1], b32 = C5_BETA[3, 2], b33 = C5_BETA[3, 3],
            b34 = C5_BETA[3, 4],
            rho1 = C5_RHO[1], rho2 = C5_RHO[2], rho3 = C5_RHO[3]),
  # plot layout and color per case
  mfrow = list(case1 = c(1, 1), case2 = c(1, 4), case3 = c(2, 4),
               case4 = c(2, 3), case5 = c(3, 6)),
  col   = list(case1 = "lightblue", case2 = "lightgreen", case3 = "lightyellow",
               case4 = "lightgreen", case5 = "lightyellow"))

# Columns differ across cases: bind within each case only.
case_of     <- vapply(raw, function(r) r$case[1L], character(1L))
all_by_case <- lapply(split(raw, case_of), function(l) do.call(rbind, l))
results <- list()

for (ci in seq_len(nrow(cells))) {
  cs <- cells$case[ci]; KK <- cells$K[ci]; TT <- cells$T[ci]; arm <- cells$arm[ci]
  res <- all_by_case[[cs]]
  res <- res[res$K == KK & res$T == TT & res$arm == arm, , drop = FALSE]
  cell_id  <- sprintf("%s_K%d_T%d%s", cs, KK, TT,
                      if (arm == "ns4") "_ns4" else "")
  cell_lab <- sprintf("K=%d,T=%d%s", KK, TT,
                      if (arm == "ns4") ",ns=4" else "")
  ok <- res$converged
  usable <- if ("label_ok" %in% names(res) && any(!is.na(res$label_ok))) {
    ok & res$label_ok & !(ok & res$bound_hit)
  } else ok
  cat(sprintf("\n=== %s | converged %d/%d | usable %d ===\n",
              cell_id, sum(ok), M, sum(usable)))

  if (sum(usable) > 0L) {
    true_vals <- TRUES[[cs]]
    hat_cols  <- paste0(names(true_vals), "_hat")
    est_mat   <- as.matrix(res[usable, hat_cols, drop = FALSE])
    colnames(est_mat) <- names(true_vals)
    smry <- make_summary(est_mat, true_vals)
    print_summary(smry)
    smry$K <- KK; smry$T <- TT; smry$arm <- arm; smry$usable <- sum(usable)
    save_csv(smry, sprintf("%s_summary.csv", cell_id))
    draw_and_save(function() boxplot_grid(est_mat, true_vals, cell_lab,
                                          TRUES$mfrow[[cs]], TRUES$col[[cs]]),
                  sprintf("%s.png", cell_id),
                  width  = 400 * TRUES$mfrow[[cs]][2] + 200,
                  height = 400 * TRUES$mfrow[[cs]][1] + 100)
  }
  results[[cell_id]] <- list(case = cs, K = KK, T = TT, arm = arm, M = M,
                             n_converged = sum(ok), n_usable = sum(usable),
                             reps = res)
}

# -- Bundle --------------------------------------------------------------------
rds_path <- file.path(results_dir, "monte_carlo_results.rds")
saveRDS(list(created = Sys.time(),
             rsdc_version = as.character(utils::packageVersion("RSDC")),
             M = M, T_grid = T_grid, T5_grid = T5_grid, K4_grid = K4_grid,
             mc_cores = MC_CORES, elapsed_min = elapsed_all / 60,
             cells = results),
        rds_path)
saved_files <- c(saved_files, rds_path)

cat("\n================================================================\n")
cat(sprintf("DONE in %.1f min. %d files written:\n", elapsed_all / 60, length(saved_files)))
for (f in saved_files) cat("  ", normalizePath(f, mustWork = FALSE), "\n", sep = "")
cat("================================================================\n")
