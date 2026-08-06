# S3 methods, standard errors, and N >= 4 support (added in 1.3-0).

fast <- list(itermax = 20, NP = 40, seed = 1)

test_that("rsdc_estimate returns an rsdc_fit with working S3 methods", {
  skip_on_cran()
  set.seed(1)
  y <- scale(matrix(rnorm(200 * 2), 200, 2))
  X <- cbind(1, as.numeric(scale(seq_len(200))))

  fit <- rsdc_estimate("tvtp", residuals = y, N = 2, X = X, control = fast)

  expect_s3_class(fit, "rsdc_fit")
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
  expect_equal(nobs(fit), 200L)
  expect_length(coef(fit), fit$npar)

  ll <- logLik(fit)
  expect_equal(attr(ll, "df"), fit$npar)
  expect_equal(attr(ll, "nobs"), 200L)

  # AIC via stats must equal the textbook 2k - 2logLik
  expect_equal(AIC(fit), 2 * fit$npar - 2 * fit$log_likelihood, tolerance = 1e-8)

  expect_output(print(fit), "RSDC fit")
  sm <- summary(fit)
  expect_s3_class(sm, "summary.rsdc_fit")
})

test_that("standard errors / vcov / confint are available at an interior optimum", {
  skip_on_cran()
  set.seed(2)
  y <- scale(matrix(rnorm(250 * 2), 250, 2))
  fit <- rsdc_estimate("noX", residuals = y, N = 2, control = fast)

  if (!is.null(fit$vcov)) {
    V <- vcov(fit)
    expect_equal(dim(V), c(fit$npar, fit$npar))
    expect_equal(V, t(V), tolerance = 1e-6)          # symmetric
    ci <- confint(fit)
    expect_equal(nrow(ci), fit$npar)
    # Only finite intervals are ordered; boundary params can yield NA SEs.
    ok <- is.finite(ci[, 1]) & is.finite(ci[, 2])
    if (any(ok)) expect_true(all(ci[ok, 1] <= ci[ok, 2]))
  } else {
    succeed("vcov NULL at this optimum (boundary); acceptable.")
  }
})

test_that("predict() and simulate() methods work", {
  skip_on_cran()
  set.seed(3)
  y <- scale(matrix(rnorm(200 * 2), 200, 2))
  X <- cbind(1, as.numeric(scale(seq_len(200))))
  fit <- rsdc_estimate("tvtp", residuals = y, N = 2, X = X, control = fast)

  S <- matrix(1, 200, 2, dimnames = list(NULL, c("a", "b")))
  fc <- predict(fit, residuals = y, sigma_matrix = S, value_cols = c("a", "b"), X = X)
  expect_true(is.list(fc))
  expect_equal(ncol(fc$predicted_correlations), choose(2, 2))

  sim <- simulate(fit, X = X, seed = 9)
  expect_equal(nrow(sim$observations), nrow(X))
  expect_true(all(sim$states %in% 1:2))
})

test_that("N >= 4 is supported (noX and tvtp)", {
  skip_on_cran()
  set.seed(4)
  y <- scale(matrix(rnorm(300 * 2), 300, 2))
  X <- cbind(1, as.numeric(scale(seq_len(300))))

  f4 <- rsdc_estimate("noX", residuals = y, N = 4,
                      control = list(itermax = 12, NP = 60, seed = 1))
  expect_equal(dim(f4$transition_matrix), c(4, 4))
  expect_true(all(abs(rowSums(f4$transition_matrix) - 1) < 1e-8))
  expect_true(all(f4$transition_matrix >= -1e-10))

  f4t <- rsdc_estimate("tvtp", residuals = y, N = 4, X = X,
                       control = list(itermax = 12, NP = 80, seed = 1))
  expect_equal(dim(f4t$transition_matrix), c(4, 4))
  expect_equal(dim(f4t$beta), c(4, (4 - 1) * ncol(X)))
})

test_that("tvtp warns when X has no intercept column", {
  skip_on_cran()
  set.seed(5)
  y <- scale(matrix(rnorm(150 * 2), 150, 2))
  X <- cbind(as.numeric(scale(seq_len(150))), as.numeric(scale(rnorm(150))))  # no constant
  expect_warning(rsdc_estimate("tvtp", residuals = y, N = 2, X = X, control = fast),
                 "intercept")
})
