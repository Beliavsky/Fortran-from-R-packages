# tests/testthat/test-random-generation.R
# Test random generation functions for all families

test_that("rgkw generates valid random samples", {
  set.seed(123)
  n <- 100
  result <- rgkw(n, alpha = 2, beta = 3, gamma = 1.5, delta = 2, lambda = 1.2)

  expect_type(result, "double")
  expect_length(result, n)
  expect_true(all(result > 0 & result < 1))
  expect_true(all(is.finite(result)))

  # Test reproducibility
  set.seed(123)
  result2 <- rgkw(n, 2, 3, 1.5, 2, 1.2)
  expect_equal(result, result2)
})

test_that("rbkw generates valid random samples", {
  set.seed(456)
  n <- 100
  result <- rbkw(n, alpha = 2, beta = 3, gamma = 1.5, delta = 2)

  expect_type(result, "double")
  expect_length(result, n)
  expect_true(all(result > 0 & result < 1))
})

test_that("rkkw generates valid random samples", {
  set.seed(789)
  n <- 100
  result <- rkkw(n, alpha = 2, beta = 3, delta = 2, lambda = 1.2)

  expect_type(result, "double")
  expect_length(result, n)
  expect_true(all(result > 0 & result < 1))
})

test_that("rekw generates valid random samples", {
  set.seed(101)
  n <- 100
  result <- rekw(n, alpha = 2, beta = 3, lambda = 1.5)

  expect_type(result, "double")
  expect_length(result, n)
  expect_true(all(result > 0 & result < 1))
})

test_that("rmc generates valid random samples", {
  set.seed(202)
  n <- 100
  result <- rmc(n, gamma = 2, delta = 3, lambda = 1.2)

  expect_type(result, "double")
  expect_length(result, n)
  expect_true(all(result > 0 & result < 1))
})

test_that("rkw generates valid random samples", {
  set.seed(303)
  n <- 100
  result <- rkw(n, alpha = 2, beta = 3)

  expect_type(result, "double")
  expect_length(result, n)
  expect_true(all(result > 0 & result < 1))
})

test_that("rbeta_ generates valid random samples", {
  set.seed(404)
  n <- 100
  result <- rbeta_(n, gamma = 2, delta = 3)

  expect_type(result, "double")
  expect_length(result, n)
  expect_true(all(result > 0 & result < 1))
})
