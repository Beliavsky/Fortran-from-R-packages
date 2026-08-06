# -----------------------------------------
# Test Script - Error for Function Inputs
# -----------------------------------------

# Required libraries
library(RPESE)

# Context of test script
context("Verify input for functions.")

# There should be an error if we want to compute the IF TS, and no returns are provided
test_that("Error for invalid function inputs", {

  # Loading data from Performance Analytics
  data(edhec, package = "PerformanceAnalytics")
  colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
                      "GM", "LS", "MA", "RV", "SS", "FoF")

  # Error for estimator function input
  expect_error(EstimatorSE(data=edhec, estimator.fun="NA", se.method="IFcor"))
  expect_error(EstimatorSE(data=edhec, estimator.fun="ES", se.method="NA"))
})
