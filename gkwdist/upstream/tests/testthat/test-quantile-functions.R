# tests/testthat/test-quantile-functions.R
# Test quantile functions for all families

test_that("qgkw works correctly", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  result <- qgkw(p, alpha = 2, beta = 3, gamma = 1.5, delta = 2, lambda = 1.2)

  expect_type(result, "double")
  expect_length(result, length(p))
  expect_true(all(result > 0 & result < 1))
  expect_true(all(is.finite(result)))

  # Test monotonicity
  expect_true(all(diff(result) > 0))

  # Test consistency with CDF
  p_test <- 0.5
  q_test <- qgkw(p_test, 2, 3, 1.5, 2, 1.2)
  p_back <- pgkw(q_test, 2, 3, 1.5, 2, 1.2)
  expect_equal(p_back, p_test, tolerance = 1e-5)
})

test_that("qbkw works correctly", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  result <- qbkw(p, alpha = 2, beta = 3, gamma = 1.5, delta = 2)

  expect_type(result, "double")
  expect_true(all(result > 0 & result < 1))
  expect_true(all(diff(result) > 0))
})

test_that("qkkw works correctly", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  result <- qkkw(p, alpha = 2, beta = 3, delta = 2, lambda = 1.2)

  expect_type(result, "double")
  expect_true(all(result > 0 & result < 1))
  expect_true(all(diff(result) > 0))
})

test_that("qekw works correctly", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  result <- qekw(p, alpha = 2, beta = 3, lambda = 1.5)

  expect_type(result, "double")
  expect_true(all(result > 0 & result < 1))
  expect_true(all(diff(result) > 0))
})

test_that("qmc works correctly", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  result <- qmc(p, gamma = 2, delta = 3, lambda = 1.2)

  expect_type(result, "double")
  expect_true(all(result > 0 & result < 1))
  expect_true(all(diff(result) > 0))
})

test_that("qkw works correctly", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  result <- qkw(p, alpha = 2, beta = 3)

  expect_type(result, "double")
  expect_true(all(result > 0 & result < 1))
  expect_true(all(diff(result) > 0))

  # Test closed-form quantile consistency
  p_test <- 0.5
  q_test <- qkw(p_test, 2, 3)
  p_back <- pkw(q_test, 2, 3)
  expect_equal(p_back, p_test, tolerance = 1e-10)
})

test_that("qbeta_ works correctly", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  result <- qbeta_(p, gamma = 2, delta = 3)

  expect_type(result, "double")
  expect_true(all(result > 0 & result < 1))
  expect_true(all(diff(result) > 0))
})
