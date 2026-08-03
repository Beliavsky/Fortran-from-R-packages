library(testthat)
library(MultiATSM)

# Simulate yield data for testing
LoadData("CM_2024")
Yields <- Yields[1:18, ]
Economies <- c("China", "Brazil", "Mexico")

# 1) Test output dimensions, type and absences of NAs
test_that("pca_weights_one_country returns correct dimensions and type", {
  J <- nrow(Yields) / length(Economies)
  for (i in 1:length(Economies)) {
    W <- pca_weights_one_country(Yields, Economies[i])
    expect_type(W, "double")
    expect_equal(dim(W), c(J, J))
    expect_false(anyNA(W))
  }
})


# 2) Test sign adjustment for level
test_that("pca_weights_one_country properly adjusts sign for level weights", {
  for (i in 1:length(Economies)) {
    W <- pca_weights_one_country(Yields, Economies[i])
    expect_true(all(W[1, ] >= 0) || all(W[1, ] <= 0))
  }
})


# 3) Test error for missing economy
test_that("pca_weights_one_country issues an error for missing economy", {
  expect_error(pca_weights_one_country(Yields, "India"))
})
