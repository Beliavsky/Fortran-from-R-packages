# tests/testthat/test-start-values.R
# Test initial values function for all families

test_that("gkwgetstartvalues works for gkw family", {
  set.seed(123)
  data <- rgkw(100, alpha = 2, beta = 3, gamma = 1.5, delta = 2, lambda = 1.2)

  result <- gkwgetstartvalues(data, family = "gkw")

  expect_type(result, "double")
  expect_length(result, 5)
  expect_true(all(result > 0))
  expect_true(all(is.finite(result)))

  # Check names
  expect_equal(names(result), c("alpha", "beta", "gamma", "delta", "lambda"))
})

test_that("gkwgetstartvalues works for bkw family", {
  set.seed(456)
  data <- rbkw(100, alpha = 2, beta = 3, gamma = 1.5, delta = 2)

  result <- gkwgetstartvalues(data, family = "bkw")

  expect_length(result, 4)
  expect_true(all(result > 0))
  expect_equal(names(result), c("alpha", "beta", "gamma", "delta"))
})

test_that("gkwgetstartvalues works for kkw family", {
  set.seed(789)
  data <- rkkw(100, alpha = 2, beta = 3, delta = 2, lambda = 1.2)

  result <- gkwgetstartvalues(data, family = "kkw")

  expect_length(result, 4)
  expect_true(all(result > 0))
  expect_equal(names(result), c("alpha", "beta", "delta", "lambda"))
})

test_that("gkwgetstartvalues works for ekw family", {
  set.seed(101)
  data <- rekw(100, alpha = 2, beta = 3, lambda = 1.5)

  result <- gkwgetstartvalues(data, family = "ekw")

  expect_length(result, 3)
  expect_true(all(result > 0))
  expect_equal(names(result), c("alpha", "beta", "lambda"))
})

test_that("gkwgetstartvalues works for mc family", {
  set.seed(202)
  data <- rmc(100, gamma = 2, delta = 3, lambda = 1.2)

  result <- gkwgetstartvalues(data, family = "mc")

  expect_length(result, 3)
  expect_true(all(result > 0))
  expect_equal(names(result), c("gamma", "delta", "lambda"))
})

test_that("gkwgetstartvalues works for kw family", {
  set.seed(303)
  data <- rkw(100, alpha = 2, beta = 3)

  result <- gkwgetstartvalues(data, family = "kw")

  expect_length(result, 2)
  expect_true(all(result > 0))
  expect_equal(names(result), c("alpha", "beta"))
})

test_that("gkwgetstartvalues works for beta family", {
  set.seed(404)
  data <- rbeta_(100, gamma = 2, delta = 3)

  result <- gkwgetstartvalues(data, family = "beta")

  expect_length(result, 2)
  expect_true(all(result > 0))
  expect_equal(names(result), c("gamma", "delta"))
})

test_that("gkwgetstartvalues handles case-insensitive family names", {
  set.seed(123)
  data <- rkw(100, alpha = 2, beta = 3)

  result1 <- gkwgetstartvalues(data, family = "kw")
  result2 <- gkwgetstartvalues(data, family = "KW")
  result3 <- gkwgetstartvalues(data, family = "Kw")

  expect_equal(result1, result2)
  expect_equal(result1, result3)
})

test_that("gkwgetstartvalues handles invalid family", {
  set.seed(123)
  data <- rkw(100, alpha = 2, beta = 3)

  expect_error(
    gkwgetstartvalues(data, family = "invalid"),
    "Invalid family"
  )
})

test_that("gkwgetstartvalues respects n_starts parameter", {
  set.seed(123)
  data <- rkw(100, alpha = 2, beta = 3)

  # Should work with different n_starts
  result1 <- gkwgetstartvalues(data, family = "kw", n_starts = 3)
  result2 <- gkwgetstartvalues(data, family = "kw", n_starts = 10)

  expect_length(result1, 2)
  expect_length(result2, 2)
  expect_true(all(result1 > 0))
  expect_true(all(result2 > 0))
})
