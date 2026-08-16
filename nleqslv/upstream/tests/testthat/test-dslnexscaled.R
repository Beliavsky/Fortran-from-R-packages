# Copyright (c) Rob Carnell 2026

test_that("scaled analytic Jacobian works for dslnex across methods and globals", {
  skip_if_not_installed("nleqslv")
  library(nleqslv)

  # Dennis–Schnabel example 6.5.1
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
  scalex <- c(2, 3)

  # First: ensure analytic Jacobian check does not fail
  z0 <- nleqslv(
    xstart, dslnex, jacdsln,
    global = "dbldog",
    control = list(btol = 0.01, delta = -1, chkjac = TRUE, scalex = scalex)
  )
  expect_false(z0$termcd == -10)

  # Broyden (default) with analytic Jacobian
  for (g in c("dbldog", "pwldog")) {
    for (delta in c(-1, -2)) {
      z <- nleqslv(
        xstart, dslnex, jacdsln,
        global = g,
        control = list(btol = 0.01, delta = delta, chkjac = TRUE, scalex = scalex)
      )
      expect_true(converged(z))
      expect_false(z$termcd == -10)
    }
  }

  # Newton with analytic Jacobian
  for (g in c("dbldog", "pwldog")) {
    for (delta in c(-1, -2)) {
      z <- nleqslv(
        xstart, dslnex, jacdsln,
        method = "Newton",
        global = g,
        control = list(btol = 0.01, delta = delta, chkjac = TRUE, scalex = scalex)
      )
      expect_true(converged(z))
      expect_false(z$termcd == -10)
    }
  }
})

