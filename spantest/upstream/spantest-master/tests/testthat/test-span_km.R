test_that("span_km returns expected structure", {
  set.seed(123)
  R1 <- matrix(rnorm(300), 100, 3)
  R2 <- matrix(rnorm(200), 100, 2)

  result <- span_km(R1, R2)

  expect_type(result, "list")
  expect_named(result, c("pval", "stat", "H0"))
  expect_type(result$stat, "double")
  expect_type(result$pval, "double")
  expect_type(result$H0, "character")

  expect_true(result$pval >= 0 && result$pval <= 1)
  expect_equal(result$H0, "delta = 0")
})

test_that("span_km handles unequal dimensions", {
  set.seed(456)
  R1 <- matrix(rnorm(500), 100, 5)
  R2 <- matrix(rnorm(100), 100, 1)

  result <- span_km(R1, R2)

  expect_true(is.numeric(result$stat))
  expect_true(length(result$pval) == 1)
})
