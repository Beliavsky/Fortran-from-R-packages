test_that("span_hk returns correct structure and types", {
  set.seed(123)
  R1 <- matrix(rnorm(300), 100, 3)
  R2 <- matrix(rnorm(200), 100, 2)

  result <- span_hk(R1, R2)

  expect_type(result, "list")
  expect_named(result, c("pval", "stat", "H0"))
  expect_type(result$stat, "double")
  expect_type(result$pval, "double")
  expect_type(result$H0, "character")

  expect_true(result$pval >= 0 && result$pval <= 1)
  expect_equal(result$H0, "alpha = 0 and delta = 0")
})

test_that("span_hk handles higher dimensions", {
  set.seed(456)
  R1 <- matrix(rnorm(1000), 100, 10)
  R2 <- matrix(rnorm(600), 100, 6)

  result <- span_hk(R1, R2)

  expect_true(is.numeric(result$stat))
  expect_true(length(result$pval) == 1)
})
