test_that("dvasicek integrates to 1", {
  f <- function(x) dvasicek(x, p = 0.03, rho = 0.10)
  result <- integrate(f, lower = 1e-10, upper = 1 - 1e-10)
  expect_equal(result$value, 1, tolerance = 1e-6)
})

test_that("pvasicek and qvasicek are inverses", {
  p <- 0.03
  rho <- 0.10
  x_vals <- c(0.01, 0.03, 0.05, 0.10, 0.20)

  for (x in x_vals) {
    prob <- pvasicek(x, p = p, rho = rho)
    x_back <- qvasicek(prob, p = p, rho = rho)
    expect_equal(x_back, x, tolerance = 1e-10)
  }
})

test_that("pvasicek CDF matches numerical integration of dvasicek", {
  p <- 0.05
  rho <- 0.15
  q_val <- 0.08

  numerical <- integrate(
    function(x) dvasicek(x, p = p, rho = rho),
    lower = 1e-10, upper = q_val
  )$value

  analytical <- pvasicek(q_val, p = p, rho = rho)
  expect_equal(analytical, numerical, tolerance = 1e-6)
})

test_that("rvasicek sample moments match theory", {
  set.seed(42)
  p <- 0.04
  rho <- 0.08
  n <- 100000
  x <- rvasicek(n, p = p, rho = rho)

  # Theoretical mean: p
  expect_equal(mean(x), p, tolerance = 0.002)

  # All values in (0,1)
  expect_true(all(x > 0 & x < 1))
})

test_that("dvasicek with log = TRUE matches log of density", {
  x <- 0.05
  d1 <- dvasicek(x, p = 0.03, rho = 0.10)
  d2 <- dvasicek(x, p = 0.03, rho = 0.10, log = TRUE)
  expect_equal(d2, log(d1))
})

test_that("pvasicek lower.tail and log.p flags work", {
  q <- 0.05
  p <- 0.03
  rho <- 0.10

  p1 <- pvasicek(q, p, rho, lower.tail = TRUE)
  p2 <- pvasicek(q, p, rho, lower.tail = FALSE)
  expect_equal(p1 + p2, 1)

  p3 <- pvasicek(q, p, rho, log.p = TRUE)
  expect_equal(p3, log(p1))
})

test_that("qvasicek lower.tail and log.p flags work", {
  p <- 0.03
  rho <- 0.10
  prob <- 0.95

  q1 <- qvasicek(prob, p, rho, lower.tail = TRUE)
  q2 <- qvasicek(1 - prob, p, rho, lower.tail = FALSE)
  expect_equal(q1, q2)

  q3 <- qvasicek(log(prob), p, rho, log.p = TRUE)
  expect_equal(q1, q3)
})

test_that("extended distribution with kappa and u works", {
  p <- 0.03
  rho <- 0.10
  kappa <- c(0.1, -0.05)
  u <- c(1.5, 0.5)

  d <- dvasicek(0.05, p, rho, kappa = kappa, u = u)
  expect_true(is.finite(d) && d > 0)

  cdf <- pvasicek(0.05, p, rho, kappa = kappa, u = u)
  expect_true(cdf > 0 && cdf < 1)

  q99 <- qvasicek(0.99, p, rho, kappa = kappa, u = u)
  expect_true(q99 > 0 && q99 < 1)

  set.seed(1)
  r <- rvasicek(100, p, rho, kappa = kappa, u = u)
  expect_length(r, 100)
  expect_true(all(r > 0 & r < 1))
})

test_that("boundary checks for p and rho", {
  expect_error(dvasicek(0.5, p = 0, rho = 0.1), "p must be in")
  expect_error(dvasicek(0.5, p = 1, rho = 0.1), "p must be in")
  expect_error(dvasicek(0.5, p = 0.5, rho = 0), "rho must be in")
  expect_error(dvasicek(0.5, p = 0.5, rho = 1), "rho must be in")
})

test_that("dvasicek returns 0 at boundaries", {
  expect_equal(dvasicek(0, p = 0.03, rho = 0.10), 0)
  expect_equal(dvasicek(1, p = 0.03, rho = 0.10), 0)
})

test_that("kappa and u length mismatch raises error", {
  expect_error(dvasicek(0.05, 0.03, 0.1, kappa = c(0.1, 0.2), u = c(1)),
               "same length")
})
