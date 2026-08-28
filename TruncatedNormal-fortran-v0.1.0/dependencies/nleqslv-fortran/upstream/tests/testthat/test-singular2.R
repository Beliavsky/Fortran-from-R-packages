# Copyright (c) Rob Carnell 2026

# Nonlinear system
f <- function(X, a, b, c1, c2, c3) {
  x <- X[1]
  y <- X[2]
  z <- X[3]

  c(
    x + y - x * y - c1,
    x + z - x * z - c2,
    a * y + b * z - c3
  )
}

# Analytic Jacobian
Jac <- function(X, a, b, c1, c2, c3) {
  x <- X[1]
  y <- X[2]
  z <- X[3]

  J <- matrix(0, nrow = 3, ncol = 3)
  J[1, 1] <- 1 - y
  J[2, 1] <- 1 - z
  J[3, 1] <- 0

  J[1, 2] <- 1 - x
  J[2, 2] <- 0
  J[3, 2] <- a

  J[1, 3] <- 0
  J[2, 3] <- 1 - x
  J[3, 3] <- b

  J
}

test_that("Newton and Broyden converge to the closed-form solution", {
  a <- 1
  b <- 1
  c1 <- 2
  c2 <- 3
  c3 <- 4

  # Closed-form solution
  x <- (a * c1 + b * c2 - c3) / (a + b - c3)
  y <- (b * c1 - b * c2 - c1 * c3 + c3) / (-a * c1 + a - b * c2 + b)
  z <- (a * (c1 - c2) + (c2 - 1) * c3) / (a * (c1 - 1) + b * (c2 - 1))
  xsol <- c(x, y, z)

  X.start <- c(1, 2, 3)

  z1 <- nleqslv(
    X.start, f, Jac,
    a = a, b = b, c1 = c1, c2 = c2, c3 = c3,
    method = "Newton",
    control = list(trace = 0, allowSingular = TRUE)
  )

  z2 <- nleqslv(
    X.start, f, Jac,
    a = a, b = b, c1 = c1, c2 = c2, c3 = c3,
    method = "Broyden",
    control = list(trace = 0, allowSingular = TRUE)
  )

  # Both solvers should match the closed-form solution
  expect_equal(z1$x, xsol, tolerance = 1e-8)
  expect_equal(z2$x, xsol, tolerance = 1e-8)

  # Function values should be near zero
  expect_true(all(abs(z1$fvec) <= 1e-8))
  expect_true(all(abs(z2$fvec) <= 1e-8))

  # Ensure no catastrophic solver failure
  expect_false(z1$termcd %in% c(-1, -2, -10))
  expect_false(z2$termcd %in% c(-1, -2, -10))
})
