library(testthat)
library(MultiATSM)

# Simulate yield data for testing
LoadData("CM_2024")
Yields <- Yields[7:18, ]
Economies <- c("Brazil", "Mexico")
N <- 3


# 1) Test correct output dimensions, labels and the absence of NAs
test_that("Spanned_Factors returns correct dimensions and labels", {
  SpaFact_TS <- Spanned_Factors(Yields, Economies, N)
  expect_type(SpaFact_TS, "double")
  expect_equal(dim(SpaFact_TS), c(length(Economies) * N, ncol(Yields)))
  expect_true(all(grepl(paste(Economies, collapse = "|"), rownames(SpaFact_TS))))
  expect_false(anyNA(SpaFact_TS))
})


# 2) Test error for N greater than number of yields
test_that("Spanned_Factors issues an error for N > number of yields", {
  expect_error(Spanned_Factors(Yields, Economies, 7))
})


# 3) Test error for invalid inputs
test_that("Spanned_Factors throws error for invalid inputs", {
  expect_error(Spanned_Factors(as.data.frame(Yields), Economies, N))
  expect_error(Spanned_Factors(Yields, character(0), N))
  expect_error(Spanned_Factors(Yields, Economies, -1))
})
