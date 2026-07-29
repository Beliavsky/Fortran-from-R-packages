simulate_fit <- function(n, p = 0.03, rho = 0.05, kappa = 0.1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  u <- rnorm(n)
  z <- rnorm(n)
  y <- pnorm((qnorm(p) + kappa * u + sqrt(rho) * z) / sqrt(1 - rho))
  vasicekfit(y ~ u, data = data.frame(y = y, u = u))
}

test_that("vcov returns symmetric PSD matrix with correct dim and names", {
  set.seed(1)
  n <- 500
  u1 <- rnorm(n); u2 <- rnorm(n); z <- rnorm(n)
  y <- pnorm((qnorm(0.03) + 0.1 * u1 - 0.05 * u2 + sqrt(0.02) * z) / sqrt(0.98))
  fit <- vasicekfit(y ~ u1 + u2, data = data.frame(y = y, u1 = u1, u2 = u2))

  V <- vcov(fit)
  expect_equal(dim(V), c(4L, 4L))
  expect_equal(rownames(V), c("p", "rho", "u1", "u2"))
  expect_equal(colnames(V), c("p", "rho", "u1", "u2"))
  expect_equal(V, t(V))
  eig <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig >= -1e-10))
})

test_that("vcov intercept-only has 2x2 matrix for (p, rho)", {
  set.seed(10)
  n <- 500
  z <- rnorm(n)
  y <- pnorm((qnorm(0.05) + sqrt(0.1) * z) / sqrt(0.9))
  fit <- vasicekfit(y ~ 1, data = data.frame(y = y))
  V <- vcov(fit)
  expect_equal(dim(V), c(2L, 2L))
  expect_equal(rownames(V), c("p", "rho"))
})

test_that("vcov rho entry matches manual delta-method calc", {
  fit <- simulate_fit(n = 500, seed = 11)
  V <- vcov(fit)
  sigma2 <- fit$sigma2
  N <- length(fit$residuals)
  m <- length(fit$kappa)
  var_sigma2 <- 2 * sigma2^2 * (N - m - 1) / N^2
  drho_dsigma2 <- 1 / (1 + sigma2)^2
  expect_equal(V["rho", "rho"], drho_dsigma2^2 * var_sigma2)
})

test_that("confint returns matrix with correct shape, names, and contains coef", {
  fit <- simulate_fit(n = 500, seed = 2)
  ci <- confint(fit)
  expect_equal(dim(ci), c(3L, 2L))
  expect_equal(rownames(ci), c("p", "rho", "u"))
  expect_equal(colnames(ci), c("2.5 %", "97.5 %"))

  cf <- coef(fit)
  expect_true(all(ci[, 1] <= cf))
  expect_true(all(cf <= ci[, 2]))
})

test_that("confint parm subsetting works by name and index", {
  fit <- simulate_fit(n = 300, seed = 3)

  ci_rho <- confint(fit, parm = "rho")
  expect_equal(dim(ci_rho), c(1L, 2L))
  expect_equal(rownames(ci_rho), "rho")

  ci_first <- confint(fit, parm = 1)
  expect_equal(rownames(ci_first), "p")

  ci_two <- confint(fit, parm = c("p", "u"))
  expect_equal(rownames(ci_two), c("p", "u"))
})

test_that("confint level controls width", {
  fit <- simulate_fit(n = 300, seed = 4)
  w95 <- confint(fit, level = 0.95)[, 2] - confint(fit, level = 0.95)[, 1]
  w99 <- confint(fit, level = 0.99)[, 2] - confint(fit, level = 0.99)[, 1]
  expect_true(all(w99 > w95))
})

test_that("SEs scale roughly as 1/sqrt(N)", {
  n1 <- 500
  se1 <- sqrt(diag(vcov(simulate_fit(n = n1,     seed = 51))))
  se2 <- sqrt(diag(vcov(simulate_fit(n = 4 * n1, seed = 52))))
  ratio <- se1 / se2
  expect_true(all(ratio > 1.5 & ratio < 2.5))
})

test_that("vcov(type = \"HAC\") returns symmetric PSD matrix with correct names", {
  skip_if_not_installed("sandwich")
  set.seed(7)
  n <- 300
  u1 <- rnorm(n); u2 <- rnorm(n); z <- rnorm(n)
  y <- pnorm((qnorm(0.03) + 0.1 * u1 - 0.05 * u2 + sqrt(0.02) * z) / sqrt(0.98))
  fit <- vasicekfit(y ~ u1 + u2, data = data.frame(y = y, u1 = u1, u2 = u2))

  V <- vcov(fit, type = "HAC")
  expect_equal(dim(V), c(4L, 4L))
  expect_equal(rownames(V), c("p", "rho", "u1", "u2"))
  expect_equal(V, t(V))
  eig <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig >= -1e-8))
})

test_that("HAC CI coverage beats iid coverage under serial correlation", {
  skip_on_cran()
  skip_if_not_installed("sandwich")

  set.seed(20260428)
  B          <- 300
  N          <- 200
  phi_u      <- 0.9
  phi_e      <- 0.5
  p_true     <- 0.03
  rho_true   <- 0.05
  kappa_true <- 0.13

  ar1 <- function(N, phi, sd = sqrt(1 - phi^2)) {
    e <- rnorm(N, sd = sd); x <- numeric(N); x[1] <- rnorm(1)
    for (i in 2:N) x[i] <- phi * x[i - 1] + e[i]
    x
  }

  cov_iid <- cov_hac <- matrix(FALSE, B, 3L,
                               dimnames = list(NULL, c("p", "rho", "u")))
  for (b in seq_len(B)) {
    u   <- ar1(N, phi_u)
    eps <- ar1(N, phi_e) * sqrt(rho_true / (1 - rho_true))
    Y   <- (qnorm(p_true) + kappa_true * u) / sqrt(1 - rho_true) + eps
    d   <- data.frame(y = pnorm(Y), u = u)
    fit <- vasicekfit(y ~ u, data = d)

    ci_i <- confint(fit, type = "iid")
    ci_h <- confint(fit, type = "HAC")
    truth <- c(p = p_true, rho = rho_true, u = kappa_true)
    for (k in seq_along(truth)) {
      cov_iid[b, k] <- ci_i[k, 1] <= truth[k] && truth[k] <= ci_i[k, 2]
      cov_hac[b, k] <- ci_h[k, 1] <= truth[k] && truth[k] <= ci_h[k, 2]
    }
  }

  rate_iid <- colMeans(cov_iid)
  rate_hac <- colMeans(cov_hac)

  # iid is badly undercovered under phi_u = 0.9
  expect_true(all(rate_iid < 0.90))
  # HAC closes most of the gap
  expect_true(all(abs(rate_hac - 0.95) < abs(rate_iid - 0.95)))
  expect_gt(rate_hac["u"] - rate_iid["u"], 0.05)
})

test_that("confint passes type through to vcov", {
  skip_if_not_installed("sandwich")
  fit <- simulate_fit(n = 300, seed = 9)
  ci_iid <- confint(fit, type = "iid")
  ci_hac <- confint(fit, type = "HAC")
  expect_equal(dim(ci_hac), dim(ci_iid))
  expect_equal(rownames(ci_hac), rownames(ci_iid))
})

test_that("vcov errors with informative message if sandwich missing", {
  skip_if(requireNamespace("sandwich", quietly = TRUE))
  fit <- simulate_fit(n = 100, seed = 10)
  expect_error(vcov(fit, type = "HAC"), "sandwich")
})

test_that("Wald CI coverage is approximately nominal under simulation", {
  skip_on_cran()
  set.seed(99)
  B <- 300
  n <- 500
  p_true <- 0.03; rho_true <- 0.05; kappa_true <- 0.1

  covered <- matrix(FALSE, B, 3L)
  for (b in seq_len(B)) {
    u <- rnorm(n); z <- rnorm(n)
    y <- pnorm((qnorm(p_true) + kappa_true * u + sqrt(rho_true) * z) /
                 sqrt(1 - rho_true))
    fit <- vasicekfit(y ~ u, data = data.frame(y = y, u = u))
    ci <- confint(fit)
    covered[b, 1] <- ci["p",   1] <= p_true     && p_true     <= ci["p",   2]
    covered[b, 2] <- ci["rho", 1] <= rho_true   && rho_true   <= ci["rho", 2]
    covered[b, 3] <- ci["u",   1] <= kappa_true && kappa_true <= ci["u",   2]
  }

  cov_rates <- colMeans(covered)
  expect_true(all(cov_rates > 0.90))
  expect_true(all(cov_rates < 0.99))
})
