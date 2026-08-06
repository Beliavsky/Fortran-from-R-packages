# Input-validation guards added in 1.4-0 (audit hardening).

test_that("rsdc_estimate rejects non-finite inputs", {
  y <- matrix(rnorm(60), 30, 2); y[5, 1] <- NA
  expect_error(rsdc_estimate("noX", residuals = y, N = 2), "NA/NaN/Inf")
  y2 <- scale(matrix(rnorm(60), 30, 2))
  X <- cbind(1, c(NA, rnorm(29)))
  expect_error(rsdc_estimate("tvtp", residuals = y2, N = 2, X = X), "NA/NaN/Inf")
})

test_that("rsdc_hamilton rejects X-without-beta and zero xi_init", {
  y <- scale(matrix(rnorm(60), 30, 2)); K <- 2; N <- 2
  rho <- rbind(0.2, 0.7)
  X <- cbind(1, as.numeric(scale(1:30)))
  expect_error(rsdc_hamilton(y, X = X, beta = NULL, rho_matrix = rho, K = K, N = N),
               "both X and beta")
  P <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, byrow = TRUE)
  expect_error(rsdc_hamilton(y, NULL, NULL, rho, K, N, P, xi_init = c(0, 0)),
               "positive sum")
})

test_that("rsdc_likelihood penalises a wrong-length parameter vector", {
  y <- scale(matrix(rnorm(60), 30, 2))
  expect_equal(rsdc_likelihood(c(0.9, 0.8), y = y, exog = NULL, K = 2, N = 2), 1e10)
})

test_that("out-of-sample estimation errors on a too-small sample", {
  skip_on_cran()
  y <- scale(matrix(rnorm(4), 2, 2))   # T = 2 -> 70/30 cut < 2
  expect_error(rsdc_estimate("const", residuals = y, out_of_sample = TRUE), "too small")
})

test_that("rsdc_estimate rejects single-column (K < 2) residuals", {
  y1 <- matrix(rnorm(30), 30, 1)
  expect_error(rsdc_estimate("noX", residuals = y1, N = 2), "at least 2 columns")
  expect_error(rsdc_estimate("const", residuals = y1), "at least 2 columns")
})

test_that("rsdc_estimate rejects X with mismatched row count", {
  y <- scale(matrix(rnorm(60), 30, 2))
  X <- cbind(1, as.numeric(scale(1:25)))   # 25 rows vs 30
  expect_error(rsdc_estimate("tvtp", residuals = y, N = 2, X = X),
               "same number of rows")
})

test_that("rsdc_hamilton rejects X with mismatched row count", {
  y <- scale(matrix(rnorm(60), 30, 2)); rho <- rbind(0.2, 0.7)
  X <- cbind(1, as.numeric(scale(1:20)))   # 20 rows vs 30
  beta <- matrix(0, 2, 2)
  expect_error(rsdc_hamilton(y, X = X, beta = beta, rho_matrix = rho, K = 2, N = 2),
               "same number of rows")
})

test_that("rsdc_hamilton rejects a non-stochastic P", {
  y <- scale(matrix(rnorm(60), 30, 2)); rho <- rbind(0.2, 0.7)
  P_bad <- matrix(c(0.9, 0.2, 0.2, 0.8), 2, byrow = TRUE)   # rows do not sum to 1
  expect_error(rsdc_hamilton(y, NULL, NULL, rho, K = 2, N = 2, P = P_bad),
               "row-stochastic")
  P_neg <- matrix(c(1.2, -0.2, 0.2, 0.8), 2, byrow = TRUE)  # negative entry
  expect_error(rsdc_hamilton(y, NULL, NULL, rho, K = 2, N = 2, P = P_neg),
               "row-stochastic")
})

test_that("rsdc_likelihood penalises K < 2", {
  y <- matrix(rnorm(30), 30, 1)
  expect_equal(rsdc_likelihood(c(0.9, 0.8, 0), y = y, exog = NULL, K = 1, N = 2), 1e10)
})
