# Copyright (c) Rob Carnell 2026

# Broyden banded function

# More, Garbow, Hillstrom: Testing Unconstrained Optimization Software
# ACM Trans. Math. Software, 7, March 1981, 17--41
# Function as shown in this paper is for optimizing not for non linear equations
# from the paper it is not clear what was used by the authors for their tests
# of non linear equation solvers
# function 31

brdban <- function(x) {
  ml <- 5
  mu <- 1
  n <- length(x)
  y <- numeric(n)

  for (k in 1:n) {
    k1 <- max(1, k - ml)
    k2 <- min(n, k + mu)

    temp <- 0
    for (j in k1:k2) {
      if (j != k) {
        temp <- temp + x[j] * (1 + x[j])
      }
    }

    y[k] <- x[k] * (2 + 5 * x[k]^2) + 1 - temp
  }
  y
}

xsol <- c(
  -0.42830, -0.47660, -0.51965, -0.55810, -0.59251,
  -0.62450, -0.62324, -0.62139, -0.62045, -0.58647
)

test_that("Known solution produces correct function values", {
  fsol <- brdban(xsol)

  expect_true(all(abs(fsol) < 1e-4))
})

test_that("nleqslv converges from xstart = -1", {
  n <- 10
  xstart <- rep(-1, n)

  znlq <- nleqslv(
    xstart, brdban,
    global = "dbldog", method = "Newton",
    control = list(trace = 0, ftol = 1e-8, xtol = 1e-8, btol = 1e-2, delta = -1)
  )

  expect_equal(znlq$termcd, 1)
  expect_equal(znlq$message, expectedMessage1)
  expect_true(all(abs(znlq$fvec) <= 1e-7))
  expect_equal(znlq$x, xsol, tolerance = 1E-5)
})

test_that("nleqslv converges from xstart = -2 with Newton method", {
  n <- 10
  xstart <- rep(-2, n)

  znlq <- nleqslv(
    xstart, brdban,
    global = "dbldog", method = "Newton",
    control = list(trace = 0, ftol = 1e-8, xtol = 1e-8, btol = 1e-2, delta = -1)
  )

  expect_equal(znlq$termcd, 1)
  expect_equal(znlq$message, expectedMessage1)
  expect_true(all(abs(znlq$fvec) <= 1e-8))
  expect_equal(znlq$x, xsol, tolerance = 1E-5)
})

test_that("nleqslv converges from xstart = -2 without specifying method", {
  n <- 10
  xstart <- rep(-2, n)

  znlq <- nleqslv(
    xstart, brdban,
    global = "dbldog",
    control = list(trace = 0, ftol = 1e-8, xtol = 1e-8, btol = 1e-2, delta = -1)
  )

  expect_equal(znlq$termcd, 1)
  expect_equal(znlq$message, expectedMessage1)
  expect_true(all(abs(znlq$fvec) <= 1e-7))
  expect_equal(znlq$x, xsol, tolerance = 1E-5)
})

test_that("Repeat convergence test for xstart = -2", {
  n <- 10
  xstart <- rep(-2, n)

  znlq <- nleqslv(
    xstart, brdban,
    global = "dbldog",
    control = list(trace = 0, ftol = 1e-8, xtol = 1e-8, btol = 1e-2, delta = -1)
  )

  expect_equal(znlq$termcd, 1)
  expect_equal(znlq$message, expectedMessage1)
  expect_true(all(abs(znlq$fvec) <= 1e-8))
  expect_equal(znlq$x, xsol, tolerance = 1E-5)
})

test_that("Broyden Banded function with testnslv", {
  n <- 10
  xstart <- rep(-1, n)
  temp <- testnslv(xstart, brdban)
  expect_true(inherits(temp, "test.nleqslv"))
  expect_true(is.data.frame(temp$out))
  expect_true(all(temp$out$termcd %in% c(1,2)))
  expect_true(all(temp$out$Iter < 25))
})

test_that("Broyden Banded with function scale", {
  temp_brdban <- function(x, fscale) brdban(x) / fscale

  n <- 10
  xstart <- -rep(1,n)

  Fscale <- rep(1,n)
  temp <- testnslv(xstart, temp_brdban, fscale=Fscale)
  expect_true(all(temp$out$Fnorm < 1E-12))
  temp <- nleqslv(xstart, temp_brdban, fscale=Fscale)

  Fscale <- c(10, 5, rep(1, n-2))
  temp2 <- testnslv(xstart, temp_brdban, fscale=Fscale)
  expect_true(all(temp2$out$Fnorm < 1E-12))
  temp2 <- nleqslv(xstart, temp_brdban, fscale=Fscale)
  expect_equal(temp$x, temp2$x, tolerance = 1E-6)

  Fscale <- c(20, 5, 2, rep(1, n-3))
  temp3 <- testnslv(xstart, temp_brdban, fscale=Fscale)
  expect_true(all(temp3$out$Fnorm < 1E-12))
  temp3 <- nleqslv(xstart, temp_brdban, fscale=Fscale)
  expect_equal(temp$x, temp3$x, tolerance = 1E-6)
})
