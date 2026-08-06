test_that("span_py returns expected structure", {
  set.seed(123)
  R1 <- matrix(rnorm(300), 100, 3)
  R2 <- matrix(rnorm(200), 100, 2)

  result <- span_py(R1, R2)

  expect_type(result, "list")
  expect_named(result, c("pval", "stat", "H0"))
  expect_type(result$pval, "double")
  expect_true(result$pval >= 0 && result$pval <= 1)
  expect_type(result$stat, "double")
  expect_type(result$H0, "character")
  expect_equal(result$H0, "alpha = 0")
})

test_that("span_py handles different dimensions correctly", {
  set.seed(42)
  R1 <- matrix(rnorm(400), 100, 4)
  R2 <- matrix(rnorm(300), 100, 3)

  result <- span_py(R1, R2)

  expect_true(is.numeric(result$stat))
  expect_true(is.numeric(result$pval))
  expect_true(length(result$stat) == 1)
  expect_true(length(result$pval) == 1)
})

test_that("span_py returns NA for a single test asset (N = 1)", {
  set.seed(1)
  R1 <- matrix(rnorm(300), 100, 3)
  R2 <- matrix(rnorm(100), 100, 1)

  result <- span_py(R1, R2)

  expect_true(is.na(result$pval))
  expect_true(is.na(result$stat))
})

