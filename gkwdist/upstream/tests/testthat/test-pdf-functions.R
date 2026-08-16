# tests/testthat/test-pdf-functions.R
# Test PDF (density) functions for all families

test_that("dgkw works correctly", {
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  result <- dgkw(x, alpha = 2, beta = 3, gamma = 1.5, delta = 2, lambda = 1.2)

  expect_type(result, "double")
  expect_length(result, length(x))
  expect_true(all(result >= 0))
  expect_true(all(is.finite(result)))

  # Test boundary behavior
  expect_equal(dgkw(0, 2, 3, 1.5, 2, 1.2), 0)
  expect_equal(dgkw(1, 2, 3, 1.5, 2, 1.2), 0)

  # Test vectorization
  expect_length(dgkw(seq(0.1, 0.9, by = 0.1), 2, 3, 1.5, 2, 1.2), 9)
})

test_that("dbkw works correctly", {
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  result <- dbkw(x, alpha = 2, beta = 3, gamma = 1.5, delta = 2)

  expect_type(result, "double")
  expect_length(result, length(x))
  expect_true(all(result >= 0))
  expect_true(all(is.finite(result)))

  expect_equal(dbkw(0, 2, 3, 1.5, 2), 0)
  expect_equal(dbkw(1, 2, 3, 1.5, 2), 0)
})

test_that("dkkw works correctly", {
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  result <- dkkw(x, alpha = 2, beta = 3, delta = 2, lambda = 1.2)

  expect_type(result, "double")
  expect_length(result, length(x))
  expect_true(all(result >= 0))
  expect_true(all(is.finite(result)))

  expect_equal(dkkw(0, 2, 3, 2, 1.2), 0)
  expect_equal(dkkw(1, 2, 3, 2, 1.2), 0)
})

test_that("dekw works correctly", {
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  result <- dekw(x, alpha = 2, beta = 3, lambda = 1.5)

  expect_type(result, "double")
  expect_length(result, length(x))
  expect_true(all(result >= 0))
  expect_true(all(is.finite(result)))

  expect_equal(dekw(0, 2, 3, 1.5), 0)
  expect_equal(dekw(1, 2, 3, 1.5), 0)
})

test_that("dmc works correctly", {
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  result <- dmc(x, gamma = 2, delta = 3, lambda = 1.2)

  expect_type(result, "double")
  expect_length(result, length(x))
  expect_true(all(result >= 0))
  expect_true(all(is.finite(result)))

  expect_equal(dmc(0, 2, 3, 1.2), 0)
  expect_equal(dmc(1, 2, 3, 1.2), 0)
})

test_that("dkw works correctly", {
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  result <- dkw(x, alpha = 2, beta = 3)

  expect_type(result, "double")
  expect_length(result, length(x))
  expect_true(all(result >= 0))
  expect_true(all(is.finite(result)))

  expect_equal(dkw(0, 2, 3), 0)
  expect_equal(dkw(1, 2, 3), 0)
})

test_that("dbeta_ works correctly", {
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  result <- dbeta_(x, gamma = 2, delta = 3)

  expect_type(result, "double")
  expect_length(result, length(x))
  expect_true(all(result >= 0))
  expect_true(all(is.finite(result)))

  expect_equal(dbeta_(0, 2, 3), 0)
  expect_equal(dbeta_(1, 2, 3), 0)
})
