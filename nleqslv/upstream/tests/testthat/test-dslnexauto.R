# Copyright (c) Rob Carnell 2026

test_that("Dennis–Schnabel example with xscalm='auto' behaves as expected", {
  library(nleqslv)

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

  # Helper: check convergence
  converged <- function(z) {
    all(abs(z$fvec) <= 1e-8)
  }

  # ------------------------------------------------------------------
  # Broyden numerical Jacobian: qline, gline
  # Both converged in your <results>
  # ------------------------------------------------------------------
  for (z in c("qline", "gline")) {
    res <- nleqslv(
      xstart, dslnex,
      global = z,
      xscalm = "auto",
      control = list(btol = 0.01)
    )
    expect_true(
      converged(res),
      info = paste("Expected convergence for global =", z, "but got:", res$message)
    )
  }

  # ------------------------------------------------------------------
  # Broyden numerical Jacobian: dbldog, pwldog
  # Only delta = -1 converged; delta = -2 stalled (per <results>)
  # ------------------------------------------------------------------
  for (z in c("dbldog", "pwldog")) {
    # delta = -1.0 → TRUE
    res1 <- nleqslv(
      xstart, dslnex,
      global = z,
      xscalm = "auto",
      control = list(btol = 0.01, delta = -1)
    )
    expect_true(
      converged(res1),
      info = paste("Expected convergence for", z, "delta=-1 but got:", res1$message)
    )

    # delta = -2.0 → FALSE (algorithm stalled)
    res2 <- nleqslv(
      xstart, dslnex,
      global = z,
      xscalm = "auto",
      control = list(btol = 0.01, delta = -2)
    )
    expect_false(
      converged(res2),
      info = paste("Expected failure (stall) for", z, "delta=-2 but got:", res2$message)
    )
  }

  # ------------------------------------------------------------------
  # Broyden analytical Jacobian
  # Same pattern: delta = -1 converged, delta = -2 stalled
  # ------------------------------------------------------------------
  for (z in c("dbldog", "pwldog")) {
    res1 <- nleqslv(
      xstart, dslnex, jacdsln,
      global = z,
      xscalm = "auto",
      control = list(btol = 0.01, delta = -1)
    )
    expect_true(
      converged(res1),
      info = paste("Expected convergence for analytic", z, "delta=-1 but got:", res1$message)
    )

    res2 <- nleqslv(
      xstart, dslnex, jacdsln,
      global = z,
      xscalm = "auto",
      control = list(btol = 0.01, delta = -2)
    )
    expect_false(
      converged(res2),
      info = paste("Expected failure (stall) for analytic", z, "delta=-2 but got:", res2$message)
    )
  }

  # ------------------------------------------------------------------
  # Newton analytical Jacobian
  # All four combinations stalled (all FALSE in <results>)
  # ------------------------------------------------------------------
  for (z in c("dbldog", "pwldog")) {
    for (delta in c(-1, -2)) {
      res <- nleqslv(
        xstart, dslnex, jacdsln,
        method = "Newton",
        global = z,
        xscalm = "auto",
        control = list(btol = 0.01, delta = delta)
      )
      expect_false(
        converged(res),
        info = paste("Expected stall for Newton analytic", z, "delta=", delta,
                     "but got:", res$message)
      )
    }
  }
})
