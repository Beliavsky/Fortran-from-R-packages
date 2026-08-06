test_that("span_grs returns correct structure and values", {
  set.seed(42)
  R1 <- matrix(rnorm(300), nrow = 100, ncol = 3)
  R2 <- matrix(rnorm(200), nrow = 100, ncol = 2)

  result <- span_grs(R1, R2)

  expect_type(result, "list")
  expect_named(result, c("pval", "stat", "H0"))
  expect_type(result$stat, "double")
  expect_type(result$pval, "double")
  expect_type(result$H0, "character")

  expect_true(result$pval >= 0 && result$pval <= 1)
  expect_equal(result$H0, "alpha = 0")
})

test_that("span_grs handles larger inputs", {
  set.seed(999)
  R1 <- matrix(rnorm(500), nrow = 100, ncol = 5)
  R2 <- matrix(rnorm(300), nrow = 100, ncol = 3)

  res <- span_grs(R1, R2)

  expect_true(is.numeric(res$stat))
  expect_length(res$pval, 1)
})
