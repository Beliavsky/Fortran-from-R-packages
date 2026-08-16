# Copyright (c) Rob Carnell 2026

# Brown almost linear function
brown <- function(x) {
  n <- length(x)
  y <- numeric(n)

  y[1:(n - 1)] <- x[1:(n - 1)] + sum(x) - (n + 1)
  y[n] <- prod(x) - 1

  y
}

# Analytic Jacobian
brownjac <- function(x) {
  n <- length(x)
  J <- matrix(1, nrow = n, ncol = n)
  diag(J) <- 2
  xprod <- prod(x)
  J[n, ] <- xprod / x
  J
}

test_that("Brown almost linear system converges for n = 50 and n = 100", {
  for (n in c(50, 100)) {
    xstart <- rep(0.5, n)

    z <- nleqslv(
      xstart,
      brown,
      jac = brownjac,
      method = "Newton",
      control = list(
        trace = 0,
        ftol = 1e-10,
        delta = "cauchy",
        allowSingular = TRUE
      )
    )

    # Original output always showed:
    #   "Function criterion near zero"
    #   TRUE
    #
    # So we assert the same condition:
    expect_true(all(abs(z$fvec) <= 1e-8))

    # And ensure solver did not fail catastrophically
    expect_equal(z$termcd, 1)
  }
})
