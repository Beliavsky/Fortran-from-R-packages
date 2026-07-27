test_that("yc_nelson_siegel fits a curve", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  expect_s3_class(fit, "yc_curve")
  expect_equal(fit$method, "nelson_siegel")
  expect_true(all(c("beta0", "beta1", "beta2", "tau") %in% names(fit$params)))
  expect_length(fit$fitted, length(d$maturities))
  expect_length(fit$residuals, length(d$maturities))
})

test_that("yc_nelson_siegel achieves good fit", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  rmse <- sqrt(mean(fit$residuals^2))
  # RMSE should be under 50 bps (0.005) for a reasonable curve

  expect_true(rmse < 0.005)
})

test_that("yc_nelson_siegel known values", {
  # Given known parameters, verify the model reproduces them
  beta0 <- 0.06
  beta1 <- -0.02
  beta2 <- 0.01
  tau <- 2

  m <- c(0.5, 1, 2, 5, 10, 30)
  # Generate synthetic data from known NS parameters
  synthetic_rates <- beta0 + beta1 * ((1 - exp(-m / tau)) / (m / tau)) +
    beta2 * ((1 - exp(-m / tau)) / (m / tau) - exp(-m / tau))

  fit <- yc_nelson_siegel(m, synthetic_rates, tau_init = 2)

  # Should recover parameters closely
  expect_equal(fit$params$beta0, beta0, tolerance = 0.001)
  expect_equal(fit$params$beta1, beta1, tolerance = 0.001)
  expect_equal(fit$params$beta2, beta2, tolerance = 0.001)
  expect_equal(fit$params$tau, tau, tolerance = 0.01)
})

test_that("yc_nelson_siegel handles flat curve", {
  d <- make_flat_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  expect_s3_class(fit, "yc_curve")
  # beta0 should be close to the flat rate
  expect_equal(fit$params$beta0, 0.04, tolerance = 0.002)
})

test_that("yc_nelson_siegel handles inverted curve", {
  d <- make_inverted_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  expect_s3_class(fit, "yc_curve")
  # beta1 should be positive (short > long)
  expect_true(fit$params$beta1 > 0)
})

test_that("yc_nelson_siegel accepts date", {
  d <- make_simple_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates, date = as.Date("2024-03-15"))
  expect_equal(fit$date, as.Date("2024-03-15"))
})

test_that("yc_nelson_siegel handles negative rates", {
  # EUR/JPY-style negative rates (2015-2022 era)
  m <- c(0.5, 1, 2, 5, 10, 30)
  rates <- c(-0.005, -0.004, -0.003, -0.001, 0.002, 0.005)
  fit <- yc_nelson_siegel(m, rates)
  expect_s3_class(fit, "yc_curve")
  rmse <- sqrt(mean(fit$residuals^2))
  expect_true(rmse < 0.005)
})

test_that("yc_nelson_siegel handles near-zero rates", {
  m <- c(1, 2, 5, 10, 30)
  rates <- c(0.001, 0.0012, 0.0015, 0.002, 0.003)
  fit <- yc_nelson_siegel(m, rates)
  expect_s3_class(fit, "yc_curve")
  rmse <- sqrt(mean(fit$residuals^2))
  expect_true(rmse < 0.005)
})

test_that("yc_nelson_siegel weighted fitting", {
  d <- make_us_yields()
  # Weight the 10Y and 30Y heavily
  w <- rep(1, length(d$maturities))
  w[d$maturities %in% c(10, 30)] <- 10
  fit_w <- yc_nelson_siegel(d$maturities, d$rates, weights = w)
  fit_uw <- yc_nelson_siegel(d$maturities, d$rates)
  # Weighted fit should have smaller residuals at 10Y and 30Y
  idx <- which(d$maturities %in% c(10, 30))
  expect_true(sum(fit_w$residuals[idx]^2) <= sum(fit_uw$residuals[idx]^2) + 1e-6)
})

test_that("yc_nelson_siegel warns with fewer than 4 points", {
  expect_warning(
    yc_nelson_siegel(c(1, 2, 5), c(0.04, 0.042, 0.043)),
    "fewer than 4"
  )
})

test_that("yc_nelson_siegel rejects bad inputs", {
  expect_error(yc_nelson_siegel("a", 0.04), "numeric")
  expect_error(yc_nelson_siegel(c(1, 2), c(0.04)), "same length")
  expect_error(yc_nelson_siegel(c(1, 2), c(0.04, 0.05), tau_init = -1), "positive")
})
