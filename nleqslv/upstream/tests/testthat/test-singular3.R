# Copyright (c) Rob Carnell 2026

# Powell cautionary example
# M.J.D. Powell, "A Hybrid Method for Nonlinear Equations",
# in Numerical Methods for Nonlinear Algebraic Equations, ed. P. Rabinowitz, 1970.

f <- function(x) {
  c(
    x[1],
    10 * x[1] / (x[1] + 0.1) + 2 * x[2]^2
  )
}

jac <- function(x) {
  J <- matrix(0, nrow = 2, ncol = 2)
  J[1, 1] <- 1
  J[1, 2] <- 0
  J[2, 1] <- 1 / (x[1] + 0.1)^2
  J[2, 2] <- 4 * x[2]
  J
}

test_that("Newton method converges for Powell example from regular and singular starts", {
  # Regular start
  z1 <- nleqslv(
    c(3, 1), f,
    method = "Newton",
    control = list(trace = 0, allowSingular = TRUE)
  )

  expect_true(all(abs(z1$fvec) <= 1e-8))
  expect_false(z1$termcd %in% c(-1, -2, -10))

  # Singular start (x2 = 0 makes Jacobian singular)
  z2 <- nleqslv(
    c(3, 0), f,
    method = "Newton",
    control = list(trace = 0, allowSingular = TRUE)
  )

  expect_true(all(abs(z2$fvec) <= 1e-8))
  expect_false(z2$termcd %in% c(-1, -2, -10))
})
