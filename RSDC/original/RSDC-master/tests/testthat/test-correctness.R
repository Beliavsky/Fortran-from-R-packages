# External correctness anchor (review Prop 2): the constant-correlation (N=1)
# Gaussian MLE with unit variances should sit very close to the Pearson sample
# correlation. This is an independent check that the likelihood/optimiser return
# sensible values, complementing the internal self-consistency tests.

test_that("const fit recovers the sample correlation (K = 2)", {
  skip_on_cran()
  set.seed(11)
  y <- scale(matrix(rnorm(400 * 2), 400, 2))   # exact unit variance after scale()

  fit <- rsdc_estimate("const", residuals = y, control = list(itermax = 80, seed = 1))
  rho_hat    <- as.numeric(fit$correlations[1, 1])
  rho_sample <- stats::cor(y)[1, 2]

  expect_true(abs(rho_hat) < 1)
  expect_equal(rho_hat, rho_sample, tolerance = 0.05)
})

test_that("const fit recovers an injected correlation (K = 2)", {
  skip_on_cran()
  set.seed(12)
  rho_true <- 0.6
  z <- matrix(rnorm(600 * 2), 600, 2)
  z[, 2] <- rho_true * z[, 1] + sqrt(1 - rho_true^2) * z[, 2]
  y <- scale(z)

  fit <- rsdc_estimate("const", residuals = y, control = list(itermax = 80, seed = 1))
  expect_equal(as.numeric(fit$correlations[1, 1]), rho_true, tolerance = 0.08)
})
