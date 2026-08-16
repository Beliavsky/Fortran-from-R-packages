# Copyright (c) Rob Carnell 2026

test_that("Cutlip reaction-rate system converges to 4 steady-state solutions", {
  # Problem 2 from Shacham (1986)
  cutlip <- function(x) {
    k1  <- 17.721
    k2  <-  3.483
    k3  <- 505.051
    kr1 <-  0.118
    kr2 <-  0.033

    r <- numeric(6)
    r[1] <- 1 - x[1] - k1 * x[1] * x[6] + kr1 * x[4]
    r[2] <- 1 - x[2] - k2 * x[2] * x[6] + kr2 * x[5]
    r[3] <- -x[3] + 2 * k3 * x[4] * x[5]
    r[4] <- k1 * x[1] * x[6] - kr1 * x[4] - k3 * x[4] * x[5]
    r[5] <- 1.5 * (k2 * x[2] * x[6] - kr2 * x[5]) - k3 * x[4] * x[5]
    r[6] <- 1 - x[4] - x[5] - x[6]
    r
  }

  # Reproducible random starts
  RNGkind(kind = "Wichmann-Hill")
  set.seed(123)

  Nrep <- 50
  xstart <- matrix(0, nrow = Nrep, ncol = 6)
  xstart[, 1] <- runif(Nrep, 0, 2)
  xstart[, 2] <- runif(Nrep, 0, 1)
  xstart[, 3] <- runif(Nrep, 0, 2)
  xstart[, 4] <- runif(Nrep, 0, 1)
  xstart[, 5] <- runif(Nrep, 0, 1)
  xstart[, 6] <- runif(Nrep, 0, 1)

  # First search
  ans <- searchZeros(
    xstart,
    cutlip,
    method = "Broyden",
    global = "dbldog"
  )

  # Expect exactly 4 distinct solutions
  expect_equal(nrow(ans$x), 4)

  # All solutions must satisfy the tolerance
  expect_true(all(ans$xfnorm <= 1e-10))

  # Restart from converged points
  zans <- searchZeros(
    ans$xstart,
    cutlip,
    method = "Broyden",
    global = "dbldog"
  )

  # Should again converge to 4 solutions
  expect_equal(length(zans$idxcvg), 4)

  # Solutions should match exactly
  expect_equal(ans$xfnorm, zans$xfnorm)
})
