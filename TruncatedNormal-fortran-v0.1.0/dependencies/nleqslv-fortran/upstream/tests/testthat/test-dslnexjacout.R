# Copyright (c) Rob Carnell 2026

test_that("Jacobian handling in nleqslv matches analytic Jacobian for dslnex", {
  skip_if_not_installed("nleqslv")
  library(nleqslv)

  # Dennis & Schnabel example 6.5.1
  dslnex <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]^2 + x[2]^2 - 2
    y[2] <- exp(x[1] - 1) + x[2]^3 - 2
    y
  }

  jacdsln <- function(x) {
    n <- length(x)
    Df <- matrix(0, n, n)
    Df[1, 1] <- 2 * x[1]
    Df[1, 2] <- 2 * x[2]
    Df[2, 1] <- exp(x[1] - 1)
    Df[2, 2] <- 3 * x[2]^2
    Df
  }

  converged <- function(z) all(abs(z$fvec) <= 1e-8)

  xstart <- c(2, 0.5)

  # 1. Internal Jacobian
  z1 <- nleqslv(xstart, dslnex, jacobian = TRUE, control = list(trace = 0))
  expect_true(converged(z1))
  expect_equal(z1$jac, jacdsln(z1$x), tolerance = 0.05)

  # 2. User-supplied Jacobian
  z2 <- nleqslv(xstart, dslnex, jacdsln, jacobian = TRUE, control = list(trace = 0))
  expect_true(converged(z2))
  expect_equal(z2$jac, jacdsln(z2$x), tolerance = 0.05)

  # 3. Newton method, internal Jacobian
  z3 <- nleqslv(xstart, dslnex, method = "Newton", jacobian = TRUE,
                control = list(trace = 0))
  expect_true(converged(z3))
  expect_equal(
    z3$jac,
    jacdsln(z3$x),
    tolerance = 1e3 * sqrt(.Machine$double.eps)
  )

  # 4. Newton method, user-supplied Jacobian
  z4 <- nleqslv(xstart, dslnex, jacdsln, method = "Newton", jacobian = TRUE,
                control = list(trace = 0))
  expect_true(converged(z4))
  expect_identical(z4$jac, jacdsln(z4$x))
})
