# Copyright (c) Rob Carnell 2026

test_that("Dennis–Schnabel example converges under multiple configurations", {

    # Problem definition
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

  xstart <- c(2, 0.5)

  ## Helper: check convergence
  expect_converged <- function(z) {
    expect_true(
      all(abs(z$fvec) <= 1e-8),
      info = paste("Message:", z$message)
    )
  }

  # ---------------------------
  # Broyden numerical Jacobian
  # ---------------------------
  for (z in c("cline", "qline", "gline")) {
    res <- nleqslv(xstart, dslnex, global = z, control = list(btol = 0.01))
    expect_converged(res)
  }

  # -------------------------------------
  # Broyden numerical Jacobian + doglegs
  # -------------------------------------
  for (z in c("dbldog", "pwldog")) {
    for (delta in c(-1.0, -2.0)) {
      res <- nleqslv(
        xstart, dslnex,
        global = z,
        control = list(btol = 0.01, delta = delta)
      )
      expect_converged(res)
    }
  }

  # -------------------------------------
  # Broyden analytical Jacobian
  # -------------------------------------
  for (z in c("dbldog", "pwldog")) {
    for (delta in c(-1.0, -2.0)) {
      res <- nleqslv(
        xstart, dslnex, jacdsln,
        global = z,
        control = list(btol = 0.01, delta = delta)
      )
      expect_converged(res)
    }
  }

  # -------------------------------------
  # Newton analytical Jacobian
  # -------------------------------------
  for (z in c("dbldog", "pwldog")) {
    for (delta in c(-1.0, -2.0)) {
      res <- nleqslv(
        xstart, dslnex, jacdsln,
        method = "Newton",
        global = z,
        control = list(btol = 0.01, delta = delta)
      )
      expect_converged(res)
    }
  }
})
