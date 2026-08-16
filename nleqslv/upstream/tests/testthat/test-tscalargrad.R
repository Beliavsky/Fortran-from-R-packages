# Copyright (c) Rob Carnell 2026

test_that("Scalar Newton method matches uniroot solutions", {
  f <- function(x) {
    2.5 * exp(-0.5 * (2 * 0.045 - x)) +
      2.5 * exp(-0.045) +
      2.5 * exp(-1.5 * x) -
      100
  }

  g1 <- function(x) {
    0.5 * 2.5 * exp(-0.5 * (2 * 0.045 - x)) -
      1.5 * 2.5 * exp(-1.5 * x)
  }

  g2 <- function(x) {
    matrix(
      0.5 * 2.5 * exp(-0.5 * (2 * 0.045 - x)) -
        1.5 * 2.5 * exp(-1.5 * x),
      nrow = 1, ncol = 1
    )
  }

  # Reference solutions from uniroot
  xu1 <- uniroot(f, c(-3, 0), tol = 1e-8)$root
  xu2 <- uniroot(f, c(6, 8),  tol = 1e-8)$root

  # nleqslv with scalar derivative
  xg1_1 <- nleqslv(-2, f, g1)$x
  xg1_2 <- nleqslv( 8, f, g1)$x

  # nleqslv with 1x1 Jacobian matrix
  xg2_1 <- nleqslv(-2, f, g2)$x
  xg2_2 <- nleqslv( 8, f, g2)$x

  # Compare to uniroot
  expect_equal(xg1_1, xu1)
  expect_equal(xg1_2, xu2)
  expect_equal(xg2_1, xu1)
  expect_equal(xg2_2, xu2)

  # Cross-check derivative forms
  expect_equal(xg1_2, xg2_2)
})
