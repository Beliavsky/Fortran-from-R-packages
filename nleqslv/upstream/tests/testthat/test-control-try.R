# Copyright (c) Rob Carnell 2026

# Dennis–Schnabel example 6.5.1 (p.149)
f <- function(x) {
  y <- numeric(2)
  y[1] <- x[1]^2 + x[2]^2 - 2
  y[2] <- exp(x[1] - 1) + x[2]^3 - 2
  y
}

test_that("control argument must be a named list", {
  expect_error(
    nleqslv(f, control = list(1e-3)),
    "named list"
  )
})

test_that("unknown names in control list are rejected", {
  expect_error(
    nleqslv(f, control = list(f = 1e-3)),
    "unknown names"
  )

  expect_error(
    nleqslv(f, control = list(f = 1e-7, b = 1e-3)),
    "unknown names"
  )
})
