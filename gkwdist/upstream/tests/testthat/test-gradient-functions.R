# tests/testthat/test-gradient-functions.R
# Test gradient (score) functions for all families

test_that("grgkw returns valid gradient", {
  set.seed(123)
  data <- rgkw(50, alpha = 2, beta = 3, gamma = 1.5, delta = 2, lambda = 1.2)
  par <- c(2, 3, 1.5, 2, 1.2)

  result <- grgkw(par, data)

  expect_type(result, "double")
  expect_length(result, 5)
  expect_true(all(is.finite(result)))
})

test_that("grbkw returns valid gradient", {
  set.seed(456)
  data <- rbkw(50, alpha = 2, beta = 3, gamma = 1.5, delta = 2)
  par <- c(2, 3, 1.5, 2)

  result <- grbkw(par, data)

  expect_type(result, "double")
  expect_length(result, 4)
  expect_true(all(is.finite(result)))
})

test_that("grkkw returns valid gradient", {
  set.seed(789)
  data <- rkkw(50, alpha = 2, beta = 3, delta = 2, lambda = 1.2)
  par <- c(2, 3, 2, 1.2)

  result <- grkkw(par, data)

  expect_type(result, "double")
  expect_length(result, 4)
  expect_true(all(is.finite(result)))
})

test_that("grekw returns valid gradient", {
  set.seed(101)
  data <- rekw(50, alpha = 2, beta = 3, lambda = 1.5)
  par <- c(2, 3, 1.5)

  result <- grekw(par, data)

  expect_type(result, "double")
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
})

test_that("grmc returns valid gradient", {
  set.seed(202)
  data <- rmc(50, gamma = 2, delta = 3, lambda = 1.2)
  par <- c(2, 3, 1.2)

  result <- grmc(par, data)

  expect_type(result, "double")
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
})

test_that("grkw returns valid gradient", {
  set.seed(303)
  data <- rkw(50, alpha = 2, beta = 3)
  par <- c(2, 3)

  result <- grkw(par, data)

  expect_type(result, "double")
  expect_length(result, 2)
  expect_true(all(is.finite(result)))
})

test_that("grbeta returns valid gradient", {
  set.seed(404)
  data <- rbeta_(50, gamma = 2, delta = 3)
  par <- c(2, 3)

  result <- grbeta(par, data)

  expect_type(result, "double")
  expect_length(result, 2)
  expect_true(all(is.finite(result)))
})
