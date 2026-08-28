# Copyright (c) Rob Carnell 2026

# High-degree polynomial system from Kearfott (1987)
hdp <- function(x) {
  f <- numeric(3)
  f[1] <- 5 * x[1]^9 - 6 * x[1]^5 * x[2]^2 + x[1] * x[2]^4 + 2 * x[1] * x[3]
  f[2] <- -2 * x[1]^6 * x[2] + 2 * x[1]^2 * x[2]^3 + 2 * x[2] * x[3]
  f[3] <- x[1]^2 + x[2]^2 - 0.265625
  f
}

test_that("searchZeros finds 12 real roots for Kearfott high-degree polynomial system", {
  set.seed(123)

  N <- 40
  xstart <- matrix(runif(3 * N, min = -1, max = 1), nrow = N, ncol = 3)

  ans <- searchZeros(
    xstart, hdp,
    method = "Broyden",
    global = "dbldog"
  )

  # Expect exactly 12 distinct real roots
  expect_equal(nrow(ans$x), 12)

  # All solutions must satisfy ||f(x)|| <= 1e-10
  expect_true(all(ans$xfnorm <= 1e-10))

  # Re-run from the discovered roots
  zans <- searchZeros(
    ans$xstart, hdp,
    method = "Broyden",
    global = "dbldog"
  )

  expect_equal(length(zans$idxcvg), 12)
})

