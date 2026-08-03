test_that("calculating effective tax rates", {
  
  test_asset_returns <- generate_test_asset_returns(9)$returns

  effective_tax_rates <- calc_effective_tax_rate(
    portfolio                = test_asset_returns,
    tax_rate_ltcg            = 0.20,
    tax_rate_ordinary_income = 0.40
  )
  if (interactive()) print(effective_tax_rates, width = Inf)
  
  expect_equal(
    round(effective_tax_rates$aftertax$effective_tax_rate, 4),
    c(0.2571, 0.2653, 0.2683, 0.2706, rep(0.4000, 5))
  )

  expect_true(
    all(
      -(effective_tax_rates$aftertax$effective_tax_rate - 1) * 
        effective_tax_rates$expected_return == 
        effective_tax_rates$aftertax$postliquidation_aftertax_expected_return
    )
  )
})
