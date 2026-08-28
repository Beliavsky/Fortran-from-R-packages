# Copyright (c) Rob Carnell 2026

test_that("Broyden banded function behaves consistently across solver options", {
  library(nleqslv)

  brdban <- function(x, ml = 5, mu = 1) {
    n <- length(x)
    y <- numeric(n)

    for (k in 1:n) {
      k1 <- max(1, k - ml)
      k2 <- min(n, k + mu)

      temp <- 0.0
      for (j in k1:k2) {
        if (j != k) {
          temp <- temp + x[j] * (1.0 + x[j])
        }
      }

      y[k] <- x[k] * (2.0 + 5.0 * x[k]^2) + 1.0 - temp
    }
    y
  }

  n <- 10
  xstart <- -rep(1, n)
  ztol <- 1000 * .Machine$double.eps

  # --- Newton method, default vs. banded Jacobian ---
  z1 <- nleqslv(xstart, brdban, method = "Newton")
  z2 <- nleqslv(xstart, brdban, method = "Newton",
                control = list(dsub = 5, dsuper = 1))

  expect_equal(z1$termcd, 1)
  expect_equal(z2$termcd, 1)
  expect_equal(z1$message, expectedMessage1)
  expect_equal(z2$message, expectedMessage1)
  expect_equal(z2$x, z1$x)
  expect_equal(z2$x, z1$x, tolerance = ztol)

  # --- Newton with ml = 2, mu = 2 ---
  z1 <- nleqslv(xstart, brdban, ml = 2, mu = 2, method = "Newton")
  z2 <- nleqslv(xstart, brdban, ml = 2, mu = 2, method = "Newton",
                control = list(dsub = 2, dsuper = 2))

  expect_equal(z1$termcd, 1)
  expect_equal(z2$termcd, 1)
  expect_equal(z1$message, expectedMessage1)
  expect_equal(z2$message, expectedMessage1)
  expect_equal(z2$x, z1$x, tolerance = ztol)

  # --- Broyden method with ml = 2, mu = 2 ---
  z3 <- nleqslv(xstart, brdban, ml = 2, mu = 2, method = "Broyden")
  z4 <- nleqslv(xstart, brdban, ml = 2, mu = 2, method = "Broyden",
                control = list(dsub = 2, dsuper = 2))

  expect_equal(z3$termcd, 1)
  expect_equal(z4$termcd, 1)
  expect_equal(z3$message, expectedMessage1)
  expect_equal(z4$message, expectedMessage1)

  expect_equal(z3$x, z1$x)
  expect_equal(z4$x, z1$x)
  expect_equal(z4$x, z3$x, tolerance = ztol)
})
