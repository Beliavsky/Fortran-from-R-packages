test_that("calculating discretionary spending with known discount rate", {

expect_equal(
    calc_discretionary_spending(
      net_worth                         = 400000,
      discount_rate                     = 0.03,
      consumption_impatience_preference = 0.02,
      smooth_consumption_preference     = 0.50,
      current_age                       = 35,
      max_age                           = 120,
      gompertz_mode                     = 90,
      gompertz_dispersion               = 10
    ),
   13167.59957,
   tolerance = 1e-5
  )
})

test_that("calculating discretionary spending with unknown discount rate", {

  test_asset_returns <- generate_test_asset_returns(9)

  test_asset_returns$returns <- 
    test_asset_returns$returns |> 
    calc_effective_tax_rate(
      tax_rate_ltcg            = 0.20,
      tax_rate_ordinary_income = 0.40
    )

    weights <- 
      test_asset_returns$returns |> 
      dplyr::select(name) |> 
      dplyr::mutate(
        human_capital = 
          c(0.0936, 0.0468, 0.0468, 0, 0.3746, 0.3746, 0, 0, 0.0635),
        liabilities = 
          c(0.1263, 0, 0, 0, 0.3368, 0.3789, 0, 0, 0.1581)
      )

  expect_equal(
    calc_discretionary_spending(
      allocations_taxable       = c(0.60, rep(0, 3), 0.40, rep(0, 4)),
      allocations_taxadvantaged = rep(0, 9),
      expected_returns    = test_asset_returns$returns$expected_return,
      standard_deviations = test_asset_returns$returns$standard_deviation,
      correlations = diag(1, 9, 9),
      income = 75000,
      nondiscretionary_consumption = 40000,
      life_insurance_premium = 0,
      effective_tax_rates = 
        test_asset_returns$returns$aftertax$effective_tax_rate,
      consumption_impatience_preference = 0.02,
      smooth_consumption_preference = 0.50,
      current_age = 25,
      max_age = 120,
      gompertz_mode = 94,
      gompertz_dispersion = 8.88,
      financial_wealth = 270500,
      liabilities = 1392064,
      liabilities_weights = weights$liabilities,
      human_capital_weights = weights$human_capital,
      human_capital = 2767689,
      risk_tolerance = 0.35
    ),
    49288.5931,
    tolerance = 1e-5
  )
})
