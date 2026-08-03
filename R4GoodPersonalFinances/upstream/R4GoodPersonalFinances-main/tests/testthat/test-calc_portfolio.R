# Test scenarios based on: 
# Idzorek, Thomas M., and Paul D. Kaplan. (2024) 
# *Lifetime Financial Advice: A Personalized Optimal Multilevel Approach*. 
# CFA Institute Research Foundation. `ISBN:978-1-952-92737-9`.
# p. 73

test_that("calculating portfolio expected return", {

  test_asset_returns <- 
    tibble::tribble(
      ~asset_class,          ~expected_return, ~standard_deviation,
      "DomesticStocks",      0.0472,           0.1588, 
      "InternationalStocks", 0.0504,           0.1718,
      "Bonds",               0.0275,           0.0562
    )
  
  if (interactive()) print(test_asset_returns)
  
  expect_error(
    calc_portfolio_expected_return(
      weights = c(0.15, 0.05, 0.99999),
      returns = test_asset_returns$expected_return
    ) 
  )
  
  expect_error(
    calc_portfolio_expected_return(
      weights = c(0.10, 0.90),
      returns = test_asset_returns$expected_return
    ) 
  )
  
  portfolio_expected_return <- calc_portfolio_expected_return(
    weights = c(0.15, 0.05, 0.80),
    returns = test_asset_returns$expected_return
  ) 
  expect_equal(
    round(portfolio_expected_return, 3),
    round(0.0319, 3)
  )
  
  portfolio_expected_return <- calc_portfolio_expected_return(
    weights = c(0.10, 0.00, 0.90),
    returns = test_asset_returns$expected_return
  ) 
  expect_equal(
    round(portfolio_expected_return, 2),
    round(0.0298, 2)
  )
  
  portfolio_expected_return <- calc_portfolio_expected_return(
    weights = c(0.30, 0.15, 0.55),
    returns = test_asset_returns$expected_return
  ) 
  expect_equal(
    round(portfolio_expected_return, 3),
    round(0.0370, 3)
  )
})
  
test_that("calculating portfolio standard deviation", {

  test_asset_returns <- 
    tibble::tribble(
      ~asset_class,          ~expected_return, ~standard_deviation,
      "DomesticStocks",      0.0472,           0.1588, 
      "InternationalStocks", 0.0504,           0.1718,
      "Bonds",               0.0275,           0.0562
    )
  if (interactive()) print(test_asset_returns)
  
  test_asset_correlations <- tibble::tribble(
    ~DomesticStocks, ~InternationalStocks, ~Bonds,
    1.00, 0.87, 0.21,
    0.87, 1.00, 0.37,
    0.21, 0.37, 1.00
  )
  test_asset_correlations <- test_asset_correlations |> as.matrix()
  
  if (interactive()) print(test_asset_correlations)
  
  portfolio_sd <- calc_portfolio_sd(
    weights             = c(0.15, 0.05, 0.80),
    standard_deviations = test_asset_returns$standard_deviation,
    correlations        = test_asset_correlations
  )
  expect_equal(
    round(portfolio_sd, 4),
    round(0.0613, 4)
  )

  portfolio_sd <- calc_portfolio_sd(
    weights             = c(0.10, 0.00, 0.90),
    standard_deviations = test_asset_returns$standard_deviation,
    correlations        = test_asset_correlations
  )
  expect_equal(
    round(portfolio_sd, 4),
    round(0.0561, 4)
  )

  portfolio_sd <- calc_portfolio_sd(
    weights             = c(0.30, 0.15, 0.55),
    standard_deviations = test_asset_returns$standard_deviation,
    correlations        = test_asset_correlations
  )
  expect_equal(
    round(portfolio_sd, 3),
    round(0.0851, 3)
  )
})

test_that("calculating portfolio parameters", {

  portfolio <- create_portfolio_template()
  expect_error(calc_portfolio_parameters(portfolio))

  portfolio$accounts$taxable <- c(10000, 30000)
  portfolio_parameters <- calc_portfolio_parameters(portfolio)

  expect_true(is.list(portfolio_parameters))
  expect_true("value" %in% names(portfolio_parameters))
  expect_true("weights" %in% names(portfolio_parameters))
  expect_true("expected_return" %in% names(portfolio_parameters))
  expect_true("standard_deviation" %in% names(portfolio_parameters))

  expect_snapshot(
    portfolio_parameters
  )
})
