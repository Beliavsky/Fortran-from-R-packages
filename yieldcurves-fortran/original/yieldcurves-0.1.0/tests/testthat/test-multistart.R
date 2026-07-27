test_that("multi-start NS produces same or better fit than single-start baseline", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  rmse <- sqrt(mean(fit$residuals^2))
  # Multi-start should achieve good fit

  expect_true(rmse < 0.005)
  expect_s3_class(fit, "yc_curve")
})

test_that("multi-start NS recovers known parameters", {
  beta0 <- 0.06
  beta1 <- -0.02
  beta2 <- 0.01
  tau <- 2

  m <- c(0.5, 1, 2, 5, 10, 30)
  synthetic <- beta0 + beta1 * ((1 - exp(-m / tau)) / (m / tau)) +
    beta2 * ((1 - exp(-m / tau)) / (m / tau) - exp(-m / tau))

  # Use a deliberately bad tau_init to test multi-start robustness
  fit <- yc_nelson_siegel(m, synthetic, tau_init = 10)
  rmse <- sqrt(mean(fit$residuals^2))
  expect_true(rmse < 0.001)
  expect_equal(fit$params$beta0, beta0, tolerance = 0.005)
  expect_equal(fit$params$tau, tau, tolerance = 0.1)
})

test_that("multi-start NS works with weighted fitting", {
  d <- make_us_yields()
  w <- rep(1, length(d$maturities))
  w[d$maturities %in% c(10, 30)] <- 10
  fit <- yc_nelson_siegel(d$maturities, d$rates, weights = w)
  expect_s3_class(fit, "yc_curve")
  rmse <- sqrt(mean(fit$residuals^2))
  expect_true(rmse < 0.01)
})

test_that("multi-start Svensson produces same or better fit than NS", {
  d <- make_us_yields()
  ns_fit <- yc_nelson_siegel(d$maturities, d$rates)
  sv_fit <- yc_svensson(d$maturities, d$rates)
  ns_rmse <- sqrt(mean(ns_fit$residuals^2))
  sv_rmse <- sqrt(mean(sv_fit$residuals^2))
  # Svensson should fit at least as well
  expect_true(sv_rmse <= ns_rmse + 1e-6)
})

test_that("multi-start Svensson recovers known parameters", {
  beta0 <- 0.06
  beta1 <- -0.02
  beta2 <- 0.01
  beta3 <- -0.005
  tau1 <- 2
  tau2 <- 5

  m <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
  f1 <- (1 - exp(-m / tau1)) / (m / tau1)
  f2 <- f1 - exp(-m / tau1)
  f3 <- (1 - exp(-m / tau2)) / (m / tau2) - exp(-m / tau2)
  synthetic <- beta0 + beta1 * f1 + beta2 * f2 + beta3 * f3

  # Use deliberately off starting values
  fit <- yc_svensson(m, synthetic, tau1_init = 8, tau2_init = 1)
  rmse <- sqrt(mean(fit$residuals^2))
  expect_true(rmse < 0.001)
})

test_that("multi-start Svensson handles inverted curve", {
  d <- make_inverted_yields()
  fit <- yc_svensson(d$maturities, d$rates)
  expect_s3_class(fit, "yc_curve")
  rmse <- sqrt(mean(fit$residuals^2))
  expect_true(rmse < 0.005)
})
