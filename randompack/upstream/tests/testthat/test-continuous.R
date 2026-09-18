test_that("unif returns numeric values in [0,1)", {
  rng <- randompack_rng("pcg64")
  x <- rng$unif(10)
  expect_type(x, "double")
  expect_length(x, 10)
  expect_true(all(is.finite(x)))
  expect_true(all(x >= 0))
  expect_true(all(x < 1))
})

test_that("normal returns numeric values with defaults", {
  rng <- randompack_rng("pcg64")
  x <- rng$normal(10)
  expect_type(x, "double")
  expect_length(x, 10)
  expect_true(all(is.finite(x)))
})

test_that("unif returns numeric values in [a,b]", {
  rng <- randompack_rng("pcg64")
  x <- rng$unif(10, -1, 2)
  expect_type(x, "double")
  expect_length(x, 10)
  expect_true(all(is.finite(x)))
  expect_true(all(x >= -1))
  expect_true(all(x <= 2))
})

test_that("other continuous distributions return numeric values", {
  rng <- randompack_rng("pcg64")
  check <- function(x, n) {
    expect_type(x, "double")
    expect_length(x, n)
    expect_true(all(is.finite(x)))
  }
  check(rng$normal(5, 0, 1), 5)
  check(rng$lognormal(5, 0, 1), 5)
  check(rng$gumbel(5, 0, 1), 5)
  check(rng$pareto(5, 1, 2), 5)
  check(rng$exp(5, 1), 5)
  check(rng$gamma(5, 2, 1), 5)
  check(rng$chi2(5, 2), 5)
  check(rng$beta(5, 2, 3), 5)
  check(rng$skew_normal(5, 0, 1, 2), 5)
  check(rng$t(5, 5), 5)
  check(rng$f(5, 5, 7), 5)
  check(rng$weibull(5, 1, 2), 5)
})
