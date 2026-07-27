test_that("yc_svensson fits a curve", {
  d <- make_us_yields()
  fit <- yc_svensson(d$maturities, d$rates)
  expect_s3_class(fit, "yc_curve")
  expect_equal(fit$method, "svensson")
  expect_true(all(c("beta0", "beta1", "beta2", "beta3", "tau1", "tau2") %in%
                    names(fit$params)))
})

test_that("yc_svensson fits at least as well as nelson_siegel", {
  d <- make_us_yields()
  ns_fit <- yc_nelson_siegel(d$maturities, d$rates)
  sv_fit <- yc_svensson(d$maturities, d$rates)

  ns_rmse <- sqrt(mean(ns_fit$residuals^2))
  sv_rmse <- sqrt(mean(sv_fit$residuals^2))

  # Svensson has more parameters, should fit at least as well
  expect_true(sv_rmse <= ns_rmse + 1e-6)
})

test_that("yc_svensson warns with fewer than 6 points", {
  expect_warning(
    yc_svensson(c(1, 2, 5, 10), c(0.04, 0.042, 0.043, 0.045)),
    "fewer than 6"
  )
})

test_that("yc_svensson known-value recovery", {
  # Generate synthetic data from known Svensson parameters
  beta0 <- 0.06
  beta1 <- -0.02
  beta2 <- 0.01
  beta3 <- -0.005
  tau1 <- 2
  tau2 <- 5

  m <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
  # Svensson formula
  f1 <- (1 - exp(-m / tau1)) / (m / tau1)
  f2 <- f1 - exp(-m / tau1)
  f3 <- (1 - exp(-m / tau2)) / (m / tau2) - exp(-m / tau2)
  synthetic_rates <- beta0 + beta1 * f1 + beta2 * f2 + beta3 * f3

  fit <- yc_svensson(m, synthetic_rates, tau1_init = 2, tau2_init = 5)
  rmse <- sqrt(mean(fit$residuals^2))
  expect_true(rmse < 0.001)
})

test_that("yc_svensson handles negative rates", {
  m <- c(0.5, 1, 2, 5, 7, 10, 20, 30)
  rates <- c(-0.005, -0.004, -0.003, -0.001, 0.001, 0.003, 0.005, 0.006)
  fit <- yc_svensson(m, rates)
  expect_s3_class(fit, "yc_curve")
  rmse <- sqrt(mean(fit$residuals^2))
  expect_true(rmse < 0.005)
})

test_that("yc_svensson rejects bad inputs", {
  expect_error(yc_svensson("a", 0.04), "numeric")
  expect_error(yc_svensson(c(1, 2), c(0.04, 0.05), tau1_init = -1), "positive")
  expect_error(yc_svensson(c(1, 2), c(0.04, 0.05), tau2_init = 0), "positive")
})
