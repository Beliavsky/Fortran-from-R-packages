# Copyright (c) Rob Carnell 2026

brdtri <- function(x) {
  n <- length(x)
  (3 - 2 * x) * x - c(0, x[-n]) - 2 * c(x[-1], 0) + 1
}

brdtrijac <- function(x) {
  n <- length(x)
  J <- diag(3 - 4 * x, n, n)
  J[row(J) == col(J) + 1] <- -1
  J[row(J) == col(J) - 1] <- -2
  J
}

n <- 10
xstart <- rep(-1, n)

# Baseline solution (no Jacobian)
z0 <- nleqslv(
  xstart, brdtri,
  method = "Newton",
  global = "dbldog"
)

# --- Tests -------------------------------------------------------------------

test_that("Baseline Newton/dbldog solve converges", {
  expect_equal(z0$termcd, 1)
  expect_equal(z0$message, expectedMessage1)
  expect_true(all(abs(z0$fvec) < 1e-7))
})

test_that("Analytic Jacobian gives same solution as baseline", {
  z1 <- nleqslv(
    xstart, brdtri, brdtrijac,
    method = "Newton",
    global = "dbldog",
    control = list(trace = 0)
  )

  expect_equal(z1$termcd, 1)
  expect_equal(z0$message, expectedMessage1)
  expect_equal(z1$x, z0$x, tolerance = 1e-12)
})

test_that("Jacobian checking enabled still matches baseline", {
  z2 <- nleqslv(
    xstart, brdtri, brdtrijac,
    method = "Newton",
    global = "dbldog",
    control = list(trace = 0, chkjac = TRUE)
  )

  expect_equal(z2$termcd, 1)
  expect_equal(z0$message, expectedMessage1)
  expect_equal(z2$x, z0$x, tolerance = 1e-12)
})

test_that("Banded Jacobian specification (dsub/dsuper) matches baseline", {
  z3 <- nleqslv(
    xstart, brdtri, brdtrijac,
    method = "Newton",
    global = "dbldog",
    control = list(trace = 0, dsub = 1, dsuper = 1)
  )

  expect_equal(z3$termcd, 1)
  expect_equal(z0$message, expectedMessage1)
  expect_equal(z3$x, z0$x, tolerance = 1e-12)
})

test_that("Banded Jacobian + chkjac matches baseline", {
  z4 <- nleqslv(
    xstart, brdtri, brdtrijac,
    method = "Newton",
    global = "dbldog",
    control = list(trace = 0, dsub = 1, dsuper = 1, chkjac = TRUE)
  )

  expect_equal(z4$termcd, 1)
  expect_equal(z0$message, expectedMessage1)
  expect_equal(z4$x, z0$x, tolerance = 1e-12)
})
