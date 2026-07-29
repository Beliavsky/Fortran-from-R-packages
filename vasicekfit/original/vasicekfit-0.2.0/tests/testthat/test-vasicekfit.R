test_that("vasicekfit recovers known parameters from large sample", {
  set.seed(123)
  n <- 5000
  p_true <- 0.03
  rho_true <- 0.05
  kappa_true <- c(unemp = 0.13, hpi = -0.07)

  u1 <- rnorm(n)
  u2 <- rnorm(n)
  z <- rnorm(n)
  x <- pnorm(
    (qnorm(p_true) + kappa_true[1] * u1 + kappa_true[2] * u2 +
      sqrt(rho_true) * z) / sqrt(1 - rho_true)
  )

  d <- data.frame(default_rate = x, unemp = u1, hpi = u2)
  fit <- vasicekfit(default_rate ~ unemp + hpi, data = d)

  expect_s3_class(fit, "vasicekfit")
  expect_equal(fit$p, p_true, tolerance = 0.005)
  expect_equal(fit$rho, rho_true, tolerance = 0.005)
  expect_equal(unname(fit$kappa["unemp"]), unname(kappa_true["unemp"]), tolerance = 0.01)
  expect_equal(unname(fit$kappa["hpi"]), unname(kappa_true["hpi"]), tolerance = 0.01)
})

test_that("intercept-only model recovers standard Vasicek", {
  set.seed(456)
  n <- 50000
  p_true <- 0.05
  rho_true <- 0.10
  z <- rnorm(n)
  x <- pnorm((qnorm(p_true) + sqrt(rho_true) * z) / sqrt(1 - rho_true))

  d <- data.frame(default_rate = x)
  fit <- vasicekfit(default_rate ~ 1, data = d)

  expect_equal(fit$p, p_true, tolerance = 0.02)
  expect_equal(fit$rho, rho_true, tolerance = 0.02)
  expect_length(fit$kappa, 0)
})

test_that("bias_correct adjusts sigma2 upward", {
  set.seed(789)
  n <- 50
  u <- rnorm(n)
  z <- rnorm(n)
  x <- pnorm((qnorm(0.03) + 0.1 * u + sqrt(0.05) * z) / sqrt(0.95))

  d <- data.frame(y = x, u = u)
  fit_raw <- vasicekfit(y ~ u, data = d, bias_correct = FALSE)
  fit_bc <- vasicekfit(y ~ u, data = d, bias_correct = TRUE)

  expect_gt(fit_bc$sigma2, fit_raw$sigma2)
})

test_that("response must be in (0, 1)", {
  d <- data.frame(y = c(0, 0.5, 1), u = 1:3)
  expect_error(vasicekfit(y ~ u, data = d), "\\(0, 1\\)")
})

test_that("coef returns named vector with p, rho, kappa", {
  set.seed(1)
  n <- 100
  u <- rnorm(n)
  x <- pnorm((qnorm(0.03) + 0.1 * u + sqrt(0.02) * rnorm(n)) / sqrt(0.98))
  d <- data.frame(y = x, u = u)
  fit <- vasicekfit(y ~ u, data = d)

  cc <- coef(fit)
  expect_named(cc, c("p", "rho", "u"))
  expect_length(cc, 3)
})

test_that("fitted and residuals are in probit space", {
  set.seed(2)
  n <- 100
  u <- rnorm(n)
  x <- pnorm((qnorm(0.03) + 0.1 * u + sqrt(0.02) * rnorm(n)) / sqrt(0.98))
  d <- data.frame(y = x, u = u)
  fit <- vasicekfit(y ~ u, data = d)

  expect_equal(unname(fitted(fit) + residuals(fit)), qnorm(x))
})

test_that("intercept-only fit matches vasicek::vsk_imm", {
  skip_if_not_installed("vasicek")

  set.seed(42)
  x <- vasicek::vsk_rvs(500, Rho = 0.1, P = 0.03)
  d <- data.frame(delq_rate = x)

  ref <- vasicek::vsk_imm(d$delq_rate)
  fit <- vasicekfit(delq_rate ~ 1, data = d)

  expect_equal(fit$p,   ref$P,   tolerance = 1e-12)
  expect_equal(fit$rho, ref$Rho, tolerance = 1e-12)
})

test_that("portfolio_size correction works", {
  set.seed(3)
  n <- 100
  u <- rnorm(n)
  x <- pnorm((qnorm(0.03) + 0.1 * u + sqrt(0.02) * rnorm(n)) / sqrt(0.98))
  d <- data.frame(y = x, u = u)

  fit1 <- vasicekfit(y ~ u, data = d)
  fit2 <- vasicekfit(y ~ u, data = d, portfolio_size = 1000)

  # With finite portfolio correction, rho should be slightly different
  expect_false(identical(fit1$rho, fit2$rho))
})

test_that("portfolio_size correction handles observed zero default rates", {
  set.seed(99)
  s <- 200
  n <- 300
  pcond <- pnorm((qnorm(0.05) + sqrt(0.08) * rnorm(n)) / sqrt(1 - 0.08))
  y <- rbinom(n, s, pcond) / s
  expect_true(any(y == 0)) # the scenario the correction exists for

  fit <- vasicekfit(y ~ 1, data = data.frame(y = y), portfolio_size = s)
  expect_true(fit$p > 0 && fit$p < 1)
  expect_true(fit$rho > 0 && fit$rho < 1)
  expect_lt(abs(fit$p - 0.05), 0.01)
  expect_lt(abs(fit$rho - 0.08), 0.03)
})

test_that("no-intercept formulas are rejected", {
  d <- data.frame(y = c(0.01, 0.02, 0.03), u = 1:3)
  expect_error(vasicekfit(y ~ u - 1, data = d), "intercept")
  expect_error(vasicekfit(y ~ u + 0, data = d), "intercept")
})
