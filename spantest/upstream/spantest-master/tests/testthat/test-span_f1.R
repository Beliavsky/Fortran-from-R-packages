test_that("span_f1 returns valid output structure and values", {
  set.seed(123)
  R1 <- matrix(rnorm(300), 100, 3)
  R2 <- matrix(rnorm(200), 100, 2)

  result <- span_f1(R1, R2)

  expect_type(result, "list")
  expect_named(result, c("pval", "stat", "H0"))
  expect_true(is.numeric(result$stat) && length(result$stat) == 1)
  expect_true(is.numeric(result$pval) && result$pval >= 0 && result$pval <= 1)
  expect_equal(result$H0, "alpha = 0")
})
