# Copyright (c) Rob Carnell 2026

test_that("hook and dbldog methods give consistent solutions for dslnex", {
  # Dennis & Schnabel example 6.5.1
  dslnex <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]^2 + x[2]^2 - 2
    y[2] <- exp(x[1] - 1) + x[2]^3 - 2
    y
  }

  # helper: check convergence
  converged <- function(z) all(abs(z$fvec) <= 1e-8)

  # helper: compare solutions
  x_equal <- function(z1, z2) all(abs(z1$x - z2$x) <= 1e-8)

  run_pair <- function(xstart) {
    # HOOK
    h1 <- nleqslv(xstart, dslnex, global = "hook",
                  control = list(btol = 0.01, delta = "cauchy", trace = 0))
    h2 <- nleqslv(xstart, dslnex, global = "hook",
                  control = list(btol = 0.01, delta = "newton", trace = 0))

    # DBLDOG
    d1 <- nleqslv(xstart, dslnex, global = "dbldog",
                  control = list(btol = 0.01, delta = "cauchy", trace = 0))
    d2 <- nleqslv(xstart, dslnex, global = "dbldog",
                  control = list(btol = 0.01, delta = "newton", trace = 0))

    list(h1 = h1, h2 = h2, d1 = d1, d2 = d2)
  }

  # ---- First starting point ----
  res1 <- run_pair(c(2, 0.5))

  expect_true(converged(res1$h1))
  expect_true(converged(res1$h2))
  expect_true(converged(res1$d1))
  expect_true(converged(res1$d2))

  expect_true(x_equal(res1$h1, res1$h2))  # hook: cauchy vs newton
  expect_true(x_equal(res1$d1, res1$d2))  # dbldog: cauchy vs newton
  expect_true(x_equal(res1$h1, res1$d1))  # hook vs dbldog

  # ---- Second starting point ----
  res2 <- run_pair(c(1.1, 1.1))

  expect_true(converged(res2$h1))
  expect_true(converged(res2$h2))
  expect_true(converged(res2$d1))
  expect_true(converged(res2$d2))

  expect_true(x_equal(res2$h1, res2$h2))  # hook: cauchy vs newton
  expect_true(x_equal(res2$d1, res2$d2))  # dbldog: cauchy vs newton
  expect_true(x_equal(res2$h1, res2$d1))  # hook vs dbldog
})
