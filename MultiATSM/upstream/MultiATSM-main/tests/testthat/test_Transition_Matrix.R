library(testthat)
library(MultiATSM)

# Set inputs
LoadData("CM_2024")
Economies <- c("China", "Brazil", "Mexico", "Uruguay")

# 1) Test Sample Mean output
test_that("Transition_Matrix returns correct output for Sample Mean", {
  W_mat <- Transition_Matrix("2000", "2002", Economies, "Sample Mean", DataConnectedness = TradeFlows)
  expect_type(W_mat, "double")
  expect_equal(dim(W_mat), c(length(Economies), length(Economies)))
  expect_equal(rownames(W_mat), Economies)
  expect_equal(colnames(W_mat), Economies)
})

# 2) Test Time-varying output
test_that("Transition_Matrix returns correct output for Time-varying", {
  W_list <- Transition_Matrix("2000", "2002", Economies, "Time-varying", DataConnectedness = TradeFlows)
  expect_type(W_list, "list")
  expect_true(all(sapply(W_list, function(x) all(dim(x) == c(length(Economies), length(Economies))))))
})


# 3) Test specific year output
test_that("Transition_Matrix returns correct output for specific year", {
  W_year <- Transition_Matrix("2000", "2002", Economies, "2001", DataConnectedness = TradeFlows)
  expect_type(W_year, "double")
  expect_equal(dim(W_year), c(length(Economies), length(Economies)))
  expect_equal(rownames(W_year), Economies)
  expect_equal(colnames(W_year), Economies)
})


# 4) Test error for missing data
test_that("Transition_Matrix returns NA matrix for missing data", {
  W_mat <- Transition_Matrix("1992", "2002", Economies, "Sample Mean", DataConnectedness = TradeFlows)
  expect_true(all(is.na(W_mat)))
})
