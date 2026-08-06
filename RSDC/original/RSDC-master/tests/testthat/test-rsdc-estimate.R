test_that("rsdc_estimate(method='const') returns coherent structure", {
  skip_on_cran()
  y <- toy_residuals(T = 12, K = 2)   # small to keep it quick

  fit <- rsdc_estimate(method = "const", residuals = y)

  expect_true(is.list(fit))
  expect_equal(dim(fit$transition_matrix), c(1, 1))
  expect_equal(ncol(fit$correlations), choose(ncol(y), 2))  # for K=2 => 1
  expect_equal(dim(fit$covariances), c(2, 2, 1))
  expect_true(is.finite(fit$log_likelihood))
})

test_that("rsdc_estimate(method='noX') returns coherent structure", {
  skip_on_cran()
  y <- toy_residuals(T = 10, K = 2)

  fit <- rsdc_estimate(method = "noX", residuals = y, N = 2)

  expect_true(is.list(fit))
  expect_equal(dim(fit$transition_matrix), c(2, 2))
  expect_equal(dim(fit$correlations), c(2, choose(ncol(y), 2)))  # 2 x 1
  expect_equal(dim(fit$covariances), c(2, 2, 2))
  expect_true(is.finite(fit$log_likelihood))
})

test_that("rsdc_estimate(method='tvtp') returns coherent structure", {
  skip_on_cran()
  y <- toy_residuals(T = 10, K = 2)
  # X is not a helper; define inline
  X <- cbind(1, scale(seq_len(nrow(y))))

  fit <- rsdc_estimate(method = "tvtp", residuals = y, N = 2, X = X)

  expect_true(is.list(fit))
  expect_equal(dim(fit$transition_matrix), c(2, 2))
  expect_equal(dim(fit$correlations), c(2, choose(ncol(y), 2)))  # 2 x 1
  expect_equal(dim(fit$covariances), c(2, 2, 2))
  expect_equal(dim(fit$beta), c(2, ncol(X)))
  expect_true(is.finite(fit$log_likelihood))
})

test_that("rsdc_estimate errors when tvtp is called without X", {
  y <- toy_residuals(T = 8, K = 2)
  expect_error(
    rsdc_estimate(method = "tvtp", residuals = y, N = 2, X = NULL),
    regexp = "X must be provided for method = 'tvtp'"
  )
})

test_that("rsdc_estimate noX N=3: returned P and correlations are consistent with IS log-likelihood", {
  skip_on_cran()
  set.seed(42)
  y <- toy_residuals(T = 30, K = 2)
  K <- ncol(y); N <- 3

  fit <- rsdc_estimate("noX", residuals = y, N = N)

  P <- fit$transition_matrix
  # N=3: 2 free trans params per row (first two columns), then N*C rho params row-wise
  trans_params <- c(P[1, 1], P[1, 2],
                    P[2, 1], P[2, 2],
                    P[3, 1], P[3, 2])
  rho_vec <- as.vector(t(fit$correlations))
  params  <- c(trans_params, rho_vec)

  nll <- rsdc_likelihood(params, y = y, exog = NULL, K = K, N = N)
  expect_equal(-nll, fit$log_likelihood, tolerance = 1e-6)
})


test_that("warm start outside the default box is kept feasible, not projected (M1)", {
  skip_on_cran()
  set.seed(42)
  y <- toy_residuals(T = 80, K = 2)

  # noX, N = 3: a valid parameter vector whose second-row free transition entry
  # (0.005) lies below the default 0.01 bound, as can happen after the
  # identifiability re-ordering of a former row-complement entry.
  start_nox <- c(0.90, 0.05,  0.005, 0.90,  0.30, 0.40,  # 6 free trans (N=3)
                 -0.20, 0.10, 0.60)                      # 3 rho (K=2)
  ll_start <- -rsdc_likelihood(start_nox, y = y, exog = NULL, K = 2, N = 3)
  expect_true(is.finite(ll_start))  # the start itself is a valid model
  fit <- rsdc_estimate("noX", residuals = y, N = 3,
                       control = list(start = start_nox, compute_se = FALSE))
  # From a feasible warm start, the local refit can only improve the likelihood.
  expect_gte(fit$log_likelihood, ll_start - 1e-6)

  # tvtp, N = 2: a beta outside the default [-10, 10] box (as produced by the
  # softmax re-referencing after relabeling) must not error or degrade the fit.
  X <- cbind(1, as.numeric(scale(seq_len(nrow(y)))))
  start_tvtp <- c(15, 0.5,  2, -0.5,   0.10, 0.60)  # beta (2x2), rho (2x1)
  ll_start2 <- -rsdc_likelihood(start_tvtp, y = y, exog = X, K = 2, N = 2)
  expect_true(is.finite(ll_start2))
  fit2 <- rsdc_estimate("tvtp", residuals = y, N = 2, X = X,
                        control = list(start = start_tvtp, compute_se = FALSE))
  expect_gte(fit2$log_likelihood, ll_start2 - 1e-6)
})
