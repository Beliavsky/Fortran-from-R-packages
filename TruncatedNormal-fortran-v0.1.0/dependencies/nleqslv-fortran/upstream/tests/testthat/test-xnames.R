# Copyright (c) Rob Carnell 2026

f <- function(x) {
  y <- numeric(length(x))
  y[1] <- x[1]^2 + x[2]^3
  y[2] <- x[1] + 2 * x[2] + 3
  y
}

test_that("nleqslv preserves names on x-values (all named)", {
  xstart <- c(a = 1.0, b = 0.5)

  z <- nleqslv(xstart, f, control = list(trace = 0))

  expect_identical(names(z$x), names(xstart))
})

test_that("nleqslv preserves names on x-values (partially named)", {
  xstart <- c(u = 1.0, 0.5)

  z <- nleqslv(xstart, f, control = list(trace = 0))

  expect_identical(names(z$x), names(xstart))
})
