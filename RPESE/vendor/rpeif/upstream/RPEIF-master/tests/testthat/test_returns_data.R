# -----------------------------------------
# Test Script - Error for returns data
# -----------------------------------------

# Required libraries
library(RPEIF)

# Context of test script
context("Verify input for functions.")

# There should be an error if we want to compute the IF TS, and no returns are provided
test_that("Error for invalid returns", {
  expect_error(IF(risk="mean", returns=NULL))
  expect_error(IF(risk="SD", returns=NULL))
  expect_error(IF(risk="SemiSD", returns=NULL))
  expect_error(IF(risk="LPM", returns=NULL))
  expect_error(IF(risk="ES", returns=NULL))
  expect_error(IF(risk="VaR", returns=NULL))
  expect_error(IF(risk="SR", returns=NULL))
  expect_error(IF(risk="SoR", returns=NULL))
  expect_error(IF(risk="ESratio", returns=NULL))
  expect_error(IF(risk="VaRratio", returns=NULL))
  expect_error(IF(risk="RachR", returns=NULL))
  expect_error(IF(risk="Omega", returns=NULL))
})