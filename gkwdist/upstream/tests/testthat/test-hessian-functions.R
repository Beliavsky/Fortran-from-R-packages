# tests/testthat/test-hessian-functions.R
# Test Hessian functions for all families

test_that("hsgkw returns valid Hessian matrix", {
  set.seed(123)
  data <- rgkw(50, alpha = 2, beta = 3, gamma = 1.5, delta = 2, lambda = 1.2)
  par <- c(2, 3, 1.5, 2, 1.2)

  result <- hsgkw(par, data)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(5, 5))
  expect_true(all(is.finite(result)))

  # Hessian should be symmetric
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("hsbkw returns valid Hessian matrix", {
  set.seed(456)
  data <- rbkw(50, alpha = 2, beta = 3, gamma = 1.5, delta = 2)
  par <- c(2, 3, 1.5, 2)

  result <- hsbkw(par, data)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(4, 4))
  expect_true(all(is.finite(result)))
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("hskkw returns valid Hessian matrix", {
  set.seed(789)
  data <- rkkw(50, alpha = 2, beta = 3, delta = 2, lambda = 1.2)
  par <- c(2, 3, 2, 1.2)

  result <- hskkw(par, data)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(4, 4))
  expect_true(all(is.finite(result)))
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("hsekw returns valid Hessian matrix", {
  set.seed(101)
  data <- rekw(50, alpha = 2, beta = 3, lambda = 1.5)
  par <- c(2, 3, 1.5)

  result <- hsekw(par, data)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 3))
  expect_true(all(is.finite(result)))
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("hsmc returns valid Hessian matrix", {
  set.seed(202)
  data <- rmc(50, gamma = 2, delta = 3, lambda = 1.2)
  par <- c(2, 3, 1.2)

  result <- hsmc(par, data)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 3))
  expect_true(all(is.finite(result)))
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("hskw returns valid Hessian matrix", {
  set.seed(303)
  data <- rkw(50, alpha = 2, beta = 3)
  par <- c(2, 3)

  result <- hskw(par, data)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2, 2))
  expect_true(all(is.finite(result)))
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("hsbeta returns valid Hessian matrix", {
  set.seed(404)
  data <- rbeta_(50, gamma = 2, delta = 3)
  par <- c(2, 3)

  result <- hsbeta(par, data)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2, 2))
  expect_true(all(is.finite(result)))
  expect_equal(result, t(result), tolerance = 1e-10)
})
