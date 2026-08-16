# tests/testthat/test-loglikelihood-functions.R
# Test log-likelihood functions for all families

test_that("llgkw works correctly", {
  set.seed(123)
  data <- rgkw(50, alpha = 2, beta = 3, gamma = 1.5, delta = 2, lambda = 1.2)
  par <- c(2, 3, 1.5, 2, 1.2)

  result <- llgkw(par, data)

  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(is.finite(result))
  expect_true(result < 0) # Log-likelihood should be negative
})

test_that("llbkw works correctly", {
  set.seed(456)
  data <- rbkw(50, alpha = 2, beta = 3, gamma = 1.5, delta = 2)
  par <- c(2, 3, 1.5, 2)

  result <- llbkw(par, data)

  expect_type(result, "double")
  expect_true(is.finite(result))
  expect_true(result < 0)
})

test_that("llkkw works correctly", {
  set.seed(789)
  data <- rkkw(50, alpha = 2, beta = 3, delta = 2, lambda = 1.2)
  par <- c(2, 3, 2, 1.2)

  result <- llkkw(par, data)

  expect_type(result, "double")
  expect_true(is.finite(result))
  expect_true(result < 0)
})

test_that("llekw works correctly", {
  set.seed(101)
  data <- rekw(50, alpha = 2, beta = 3, lambda = 1.5)
  par <- c(2, 3, 1.5)

  result <- llekw(par, data)

  expect_type(result, "double")
  expect_true(is.finite(result))
  expect_true(result < 0)
})

test_that("llmc works correctly", {
  set.seed(202)
  data <- rmc(50, gamma = 2, delta = 3, lambda = 1.2)
  par <- c(2, 3, 1.2)

  result <- llmc(par, data)

  expect_type(result, "double")
  expect_true(is.finite(result))
  expect_true(result < 0)
})

test_that("llkw works correctly", {
  set.seed(303)
  data <- rkw(50, alpha = 2, beta = 3)
  par <- c(2, 3)

  result <- llkw(par, data)

  expect_type(result, "double")
  expect_true(is.finite(result))
  expect_true(result < 0)
})

test_that("llbeta works correctly", {
  set.seed(404)
  data <- rbeta_(50, gamma = 2, delta = 3)
  par <- c(2, 3)

  result <- llbeta(par, data)

  expect_type(result, "double")
  expect_true(is.finite(result))
  expect_true(result < 0)
})
