test_that("calculate utility for parameter != 1", {

  risk_aversion  <- 1
  risk_tolerance <- 1 / risk_aversion
  if (interactive()) print(risk_tolerance)
  expect_equal(calc_utility(-1,  risk_tolerance),  0)
  expect_equal(calc_utility(0,   risk_tolerance),  0)
  expect_equal(calc_utility(1,   risk_tolerance),  0)
  expect_equal(calc_utility(0.5, risk_tolerance), -0.69314718)
  expect_equal(calc_utility(1.5, risk_tolerance),  0.405465108)
  
  risk_aversion  <- 2
  risk_tolerance <- 1 / risk_aversion
  if (interactive()) print(risk_tolerance)
  expect_equal(calc_utility(-1,  risk_tolerance),  0)
  expect_equal(calc_utility(0,   risk_tolerance),  0)
  expect_equal(calc_utility(1,   risk_tolerance),  0)
  expect_equal(calc_utility(0.5, risk_tolerance), -1)
  expect_equal(calc_utility(2,   risk_tolerance),  0.5)
  
  risk_aversion <- 3
  risk_tolerance <- 1 / risk_aversion
  if (interactive()) print(risk_tolerance)
  expect_equal(calc_utility(-1,  risk_tolerance),  0)
  expect_equal(calc_utility(0,   risk_tolerance),  0)
  expect_equal(calc_utility(1,   risk_tolerance),  0)
  expect_equal(calc_utility(0.5, risk_tolerance), -1.5)
  expect_equal(calc_utility(2,   risk_tolerance),  0.375)

  risk_tolerance <- 0
  expect_equal(calc_utility(-1,  risk_tolerance), 0)
  expect_equal(calc_utility(0,   risk_tolerance), 0)
  expect_equal(calc_utility(1,   risk_tolerance), 0)
  expect_equal(calc_utility(0.5, risk_tolerance), 0)
  expect_equal(calc_utility(2,   risk_tolerance), 0)

  risk_aversion <- 100
  risk_tolerance <- 1 / risk_aversion
  if (interactive()) print(risk_tolerance)
  expect_equal(calc_utility(-1,  risk_tolerance),  0)
  expect_equal(calc_utility(0,   risk_tolerance),  0)
  expect_equal(calc_utility(1,   risk_tolerance),  0)
  expect_equal(calc_utility(2,   risk_tolerance),  0.01010101)
})

test_that("calculate utility for parameter == 1", {

  risk_aversion  <- 1
  risk_tolerance <- 1 / risk_aversion
  if (interactive()) print(risk_tolerance)
  expect_equal(calc_utility(-1,   risk_tolerance),   0)
  expect_equal(calc_utility(0,   risk_tolerance),   0)
  expect_equal(calc_utility(1,   risk_tolerance),   0)
  expect_equal(round(calc_utility(0.5, risk_tolerance), 2), -0.69)
  expect_equal(round(calc_utility(1.5, risk_tolerance), 2),  0.41)
})

test_that("calculate utility when spending is zero", {

  expect_equal(calc_utility(0,   0),   0)
  expect_equal(calc_utility(0,   1/2), 0)
  expect_equal(calc_utility(0,   1),   0)
})

test_that("calculate utility when x is negative", {

  expect_equal(calc_utility(x = -1, 1/2), 0)
  expect_equal(calc_utility(x = -1, 0), 0)
  expect_equal(calc_utility(x = -1, 1), 0)
})

test_that("calculate utility when x is a vector with negative or zero", {

  expect_equal(
    calc_utility(
      x         = c(exp(1), 0, -1, 0, -1),
      parameter = 1/ c(1, 2, 2, 1,  1)
    ),
    c(1, 0, 0, 0,  0))
})

test_that("calculate utility on a vector", {

  risk_aversion  <- c(2, 2, 2)
  risk_tolerance <- 1 / risk_aversion
  spending       <- c(1, 0.5, 2)
  expect_equal(
    calc_utility(spending, risk_tolerance),
    c(0, -1, 0.5))
})

test_that("calculating inverse utility", {

  consumption <- 1
  
  parameter <- 0
  expect_equal(
    calc_inverse_utility(calc_utility(consumption, parameter), parameter),
    consumption
  )

  parameter <- 0.0001
  expect_equal(
    calc_inverse_utility(calc_utility(consumption, parameter), parameter),
    consumption
  )

  parameter <- 0.1
  expect_equal(
    calc_inverse_utility(calc_utility(consumption, parameter), parameter),
    consumption
  )

  parameter <- 0.2
  expect_equal(
    calc_inverse_utility(calc_utility(consumption, parameter), parameter),
    consumption
  )

  parameter <- 0.9999
  expect_equal(
    calc_inverse_utility(calc_utility(consumption, parameter), parameter),
    consumption
  )
  
  parameter <- 1
  expect_equal(
    calc_inverse_utility(calc_utility(consumption, parameter), parameter),
    consumption
  )
})

test_that("calculating inverse utility for a vector", {

  consumption <- c(1, 2)
  if (interactive()) print(consumption)
  
  parameter <- 0.1
  expect_equal(
    calc_inverse_utility(calc_utility(consumption, parameter), parameter),
    consumption
  )

  parameter <- 0.2
  expect_equal(
    calc_inverse_utility(calc_utility(consumption, parameter), parameter),
    consumption
  )

  parameter <- 0.9999
  expect_equal(
    calc_inverse_utility(calc_utility(consumption, parameter), parameter),
    consumption
  )
  
  parameter <- 1
  expect_equal(
    calc_inverse_utility(calc_utility(consumption, parameter), parameter),
    consumption
  )
})
