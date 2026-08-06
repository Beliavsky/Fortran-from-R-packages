# -----------------------------------------
# Test Script - Error for Input Data
# -----------------------------------------

# Required libraries
library(RPEGLMEN)

# Context of test script
context("Verify input for functions.")

# There should be an error if we want to compute the IF TS, and no returns are provided
test_that("Error for invalid input data", {
  expect_error(glmnet_exp())
  expect_error(glmnet_exp(A=rnorm(10), b=rnorm(10)))
  expect_error(glmnet_exp(A=matrix(rnorm(10), nrow=2), b=matrix(rnorm(10), nrow=2)))
  expect_error(fit.glmGammaNet())
  expect_error(fit.glmGammaNet(A=rnorm(10), b=rnorm(10)))
  expect_error(fit.glmGammaNet(A=matrix(rnorm(10), nrow=2), b=matrix(rnorm(10), nrow=2)))
})