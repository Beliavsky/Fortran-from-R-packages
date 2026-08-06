# Regression tests for the guards added after the 1.7-0 audits. Each of these
# reproduces a defect that shipped silently at some point during development.

test_that("a converged polish step is not reported as optimizer failure 52", {
  skip_on_cran()
  # The global stage already returns an L-BFGS-B-refined point, so the polish
  # step starts at an optimum and L-BFGS-B reports 52 with no improvement.
  # That is benign and must not surface as a non-zero convergence code.
  sim <- rsdc_simulate(n = 400, X = matrix(1, 400, 1),
                       beta = matrix(stats::qlogis(0.9), 2, 1),
                       mu = matrix(0, 2, 2),
                       sigma = array(c(1, 0.1, 0.1, 1,
                                       1, 0.8, 0.8, 1), c(2, 2, 2)),
                       N = 2, seed = 2)
  fit <- rsdc_estimate("noX", residuals = sim$observations, N = 2,
                       control = list(itermax = 40, seed = 1))
  expect_identical(as.integer(fit$convergence), 0L)
  expect_true(is.finite(as.numeric(logLik(fit))))
})

test_that("N must be a whole number", {
  set.seed(1)
  y <- matrix(rnorm(200), ncol = 2)
  expect_error(rsdc_estimate("noX", residuals = y, N = 2.7), "whole number")
  expect_error(rsdc_estimate("tvtp", residuals = y, N = 2.5,
                             X = cbind(1, seq_len(100))), "whole number")
})

test_that("standard errors survive a legitimate high-curvature fit", {
  skip_on_cran()
  # |Hessian| is O(T) and O((1 - rho^2)^-2), so a long, highly correlated
  # sample legitimately exceeds any fixed magnitude cutoff. Such a fit must
  # keep its standard errors; only a penalty-contaminated one loses them.
  set.seed(4)
  mvn <- function(n, rho) {
    L <- chol(matrix(c(1, rho, rho, 1), 2))
    matrix(rnorm(2 * n), n, 2) %*% L
  }
  y <- rbind(mvn(5000, 0.3), mvn(5000, 0.995))
  fn <- function(p) rsdc_likelihood(p, y = y, exog = NULL, K = 2, N = 2)
  par <- c(0.98, 0.98, 0.3, 0.995)
  expect_gt(max(abs(stats::optimHess(par, fn))), 1e7)   # genuinely large
  V <- RSDC:::.rsdc_vcov(par, fn)
  expect_false(is.null(V))                              # but not refused

  # A correlation pinned against its bound is still refused, with a warning.
  expect_warning(V2 <- RSDC:::.rsdc_vcov(c(0.9, 0.9, 0.3, 0.9995), fn),
                 "finite-difference|feasible")
  expect_null(V2)
})

test_that("portfolio inputs must form a positive-definite correlation matrix", {
  skip_on_cran()
  set.seed(1)
  Tn <- 8; K <- 3; nm <- c("a", "b", "c")
  sig <- matrix(1, Tn, K, dimnames = list(NULL, nm))
  ret <- matrix(rnorm(Tn * K), Tn, K, dimnames = list(NULL, nm))
  # Pairwise legal, jointly indefinite: without the check, rsdc_maxdiv()
  # reported a diversification ratio above the sqrt(K) bound.
  bad <- matrix(rep(c(0.9, 0.9, -0.9), each = Tn), Tn, 3)
  expect_error(rsdc_minvar(sig, nm, bad, ret), "positive-definite")
  expect_error(rsdc_maxdiv(sig, nm, bad, ret), "positive-definite")

  ok <- matrix(rep(c(0.3, 0.4, 0.2), each = Tn), Tn, 3)
  md <- rsdc_maxdiv(sig, nm, ok, ret, lag = TRUE)
  expect_lte(max(md$diversification_ratios), sqrt(K) + 1e-8)
  expect_identical(md$n_fallback, 0L)
})

test_that("optimizer fallbacks are counted and reported", {
  skip_on_cran()
  set.seed(1)
  Tn <- 6; K <- 2; nm <- c("a", "b")
  sig <- matrix(1, Tn, K, dimnames = list(NULL, nm))
  ret <- matrix(rnorm(Tn * K), Tn, K, dimnames = list(NULL, nm))
  pc  <- matrix(0.5, Tn, 1)
  mv <- rsdc_minvar(sig, nm, pc, ret, lag = TRUE)
  expect_identical(mv$n_fallback, 0L)
  expect_true("n_fallback" %in% names(mv))
})

test_that("forecast paths reject malformed covariance inputs and covariates", {
  skip_on_cran()
  set.seed(1)
  Tn <- 40; K <- 2
  y <- matrix(rnorm(Tn * K), Tn, K)
  fit <- rsdc_estimate("noX", residuals = y, N = 2,
                       control = list(itermax = 25, seed = 1))
  sig <- matrix(1, Tn, K, dimnames = list(NULL, c("a", "b")))
  sig[1, 1] <- -1                                   # impossible volatility
  expect_error(
    rsdc_forecast("noX", N = 2, residuals = y, final_params = fit,
                  sigma_matrix = sig, value_cols = c("a", "b")),
    "positive standard deviations")

  X <- cbind(1, as.numeric(scale(seq_len(Tn))))
  ft <- rsdc_estimate("tvtp", residuals = y, N = 2, X = X,
                      control = list(itermax = 25, seed = 1))
  expect_error(
    rsdc_forecast_ahead(ft, horizon = 2,
                        X_future = matrix(c(1, NA, 1, 0), 2, 2, byrow = TRUE)),
    "NA/NaN/Inf")
})
