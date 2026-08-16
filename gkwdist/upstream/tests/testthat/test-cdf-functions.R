# tests/testthat/test-cdf-functions.R
# Test CDF (cumulative distribution) functions for all families

test_that("pgkw works correctly", {
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  result <- pgkw(x, alpha = 2, beta = 3, gamma = 1.5, delta = 2, lambda = 1.2)

  expect_type(result, "double")
  expect_length(result, length(x))
  expect_true(all(result >= 0 & result <= 1))
  expect_true(all(is.finite(result)))

  # Test monotonicity
  expect_true(all(diff(result) >= 0))

  # Test boundaries
  expect_equal(pgkw(0, 2, 3, 1.5, 2, 1.2), 0)
  expect_equal(pgkw(1, 2, 3, 1.5, 2, 1.2), 1)
})

test_that("pbkw works correctly", {
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  result <- pbkw(x, alpha = 2, beta = 3, gamma = 1.5, delta = 2)

  expect_type(result, "double")
  expect_true(all(result >= 0 & result <= 1))
  expect_true(all(diff(result) >= 0))
})

test_that("pkkw works correctly", {
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  result <- pkkw(x, alpha = 2, beta = 3, delta = 2, lambda = 1.2)

  expect_type(result, "double")
  expect_true(all(result >= 0 & result <= 1))
  expect_true(all(diff(result) >= 0))
})

test_that("pekw works correctly", {
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  result <- pekw(x, alpha = 2, beta = 3, lambda = 1.5)

  expect_type(result, "double")
  expect_true(all(result >= 0 & result <= 1))
  expect_true(all(diff(result) >= 0))
})

test_that("pmc works correctly", {
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  result <- pmc(x, gamma = 2, delta = 3, lambda = 1.2)

  expect_type(result, "double")
  expect_true(all(result >= 0 & result <= 1))
  expect_true(all(diff(result) >= 0))
})

test_that("pkw works correctly", {
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  result <- pkw(x, alpha = 2, beta = 3)

  expect_type(result, "double")
  expect_true(all(result >= 0 & result <= 1))
  expect_true(all(diff(result) >= 0))
})

test_that("pbeta_ works correctly", {
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  result <- pbeta_(x, gamma = 2, delta = 3)

  expect_type(result, "double")
  expect_true(all(result >= 0 & result <= 1))
  expect_true(all(diff(result) >= 0))
})
