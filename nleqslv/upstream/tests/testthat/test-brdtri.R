# Copyright (c) Rob Carnell 2026

# More, Garbow, Hillstrom: Testing Unconstrained Optimization Software
# ACM Trans. Math. Software, 7, March 1981, 17--41
# Function as shown in this paper is for optimizing not for non linear equations
# from the paper it is not clear what was used by the authors for their tests
# of non linear equation solvers
# function 30

# Broyden tridiagonal function
brdtri <- function(x) {
  n <- length(x)
  y <- numeric(n)

  y[1] <- (3 - 2 * x[1]) * x[1] - 2 * x[2] + 1
  y[n] <- (3 - 2 * x[n]) * x[n] - x[n - 1] + 1

  k <- 2:(n - 1)
  y[k] <- (3 - 2 * x[k]) * x[k] - x[k - 1] - 2 * x[k + 1] + 1

  y
}

test_that("Newton method: structured vs unstructured Jacobian give same solution", {
  n <- 100
  xstart <- rep(-1, n)
  ztol <- 1000 * .Machine$double.eps

  z1 <- nleqslv(xstart, brdtri, method = "Newton")
  z2 <- nleqslv(
    xstart, brdtri,
    method = "Newton",
    control = list(dsub = 1, dsuper = 1)
  )

  expect_equal(z1$termcd, 1)
  expect_equal(z2$termcd, 1)

  expect_equal(z1$njcnt, 4)
  expect_equal(z2$njcnt, 4)

  expect_equal(z1$nfcnt, 4)
  expect_equal(z2$nfcnt, 4)

  expect_equal(z1$message, expectedMessage1)
  expect_equal(z2$message, expectedMessage1)

  expect_equal(z2$x, z1$x, tolerance = ztol)
})

test_that("Newton method repeat run: structured vs unstructured Jacobian", {
  n <- 100
  xstart <- rep(-1, n)
  ztol <- 1000 * .Machine$double.eps

  z1 <- nleqslv(xstart, brdtri, method = "Newton")
  z2 <- nleqslv(
    xstart, brdtri,
    method = "Newton",
    control = list(dsub = 1, dsuper = 1)
  )

  expect_equal(z1$termcd, 1)
  expect_equal(z2$termcd, 1)

  expect_equal(z1$njcnt, 4)
  expect_equal(z2$njcnt, 4)

  expect_equal(z1$nfcnt, 4)
  expect_equal(z2$nfcnt, 4)

  expect_equal(z1$message, expectedMessage1)
  expect_equal(z2$message, expectedMessage1)

  expect_equal(z2$x, z1$x, tolerance = ztol)
})

test_that("Broyden method: structured vs unstructured Jacobian give same solution", {
  n <- 100
  xstart <- rep(-1, n)
  ztol <- 1000 * .Machine$double.eps

  z1 <- nleqslv(xstart, brdtri, method = "Newton")  # reference solution
  z3 <- nleqslv(xstart, brdtri, method = "Broyden")
  z4 <- nleqslv(
    xstart, brdtri,
    method = "Broyden",
    control = list(dsub = 1, dsuper = 1)
  )

  expect_equal(z3$termcd, 2)
  expect_equal(z4$termcd, 2)

  expect_equal(z3$njcnt, 1)
  expect_equal(z4$njcnt, 1)

  expect_equal(z3$nfcnt, 10)
  expect_equal(z4$nfcnt, 10)

  expect_equal(z3$message, expectedMessage2)
  expect_equal(z4$message, expectedMessage2)

  # Compare Broyden solutions to Newton reference
  expect_equal(z3$x, z1$x, tolerance = sqrt(.Machine$double.eps))
  expect_equal(z4$x, z1$x, tolerance = sqrt(.Machine$double.eps))

  # Compare structured vs unstructured Broyden
  expect_equal(z4$x, z3$x, tolerance = ztol)
})
