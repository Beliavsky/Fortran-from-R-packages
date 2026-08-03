test_that("calculating expected utility", {

  expect_equal(
    calc_expected_utility(
      expected_return = 0.04,
      variance        = 0.15,
      risk_tolerance  = 0.50
    ),
    -1.0948879
  )

  expect_equal(
    calc_expected_utility(
      expected_return = 0.05,
      variance        = 0.15,
      risk_tolerance  = 0.50
    ),
    -1.08195659
  )
})


test_that("calculating optimal MVO allocations for 2 assets", {

  skip_on_cran()

  n <- 2
  test_asset_returns <- generate_test_asset_returns(n)$returns
  if (interactive()) print(test_asset_returns)
  test_asset_correlations <- generate_test_asset_returns(n)$correlations
  if (interactive()) print(test_asset_correlations)

  optimal_allocations <- 
    calc_optimal_portfolio(
      risk_tolerance      = 0.50,
      expected_returns    = test_asset_returns$expected_return,
      standard_deviations = test_asset_returns$standard_deviation,
      correlations        = test_asset_correlations,
      asset_names         = test_asset_returns$name
   )$allocations
  if (interactive()) print(optimal_allocations)

  expect_equal(
    round(optimal_allocations$total, 4),
    c(0.58480, 0.41520)
  )
  expect_equal(sum(optimal_allocations$total), 1)
})

test_that("calculating optimal MVO allocations for 3 assets", {
  
  skip_on_cran()

  n <- 3
  test_asset_returns <- generate_test_asset_returns(n)$returns
  if (interactive()) print(test_asset_returns)
  test_asset_correlations <- generate_test_asset_returns(n)$correlations
  if (interactive()) print(test_asset_correlations)

  optimal_allocations <- 
    calc_optimal_portfolio(
    risk_tolerance      = 0.35,
    expected_returns    = test_asset_returns$expected_return,
    standard_deviations = test_asset_returns$standard_deviation,
    correlations        = test_asset_correlations
  )$allocations 
  if (interactive()) print(optimal_allocations)

  expect_equal(
    round(optimal_allocations$total, 3),
    c(0.273, 0.089, 0.638)
  )
  expect_equal(sum(optimal_allocations$total), 1)
})

test_that("calculating optimal MVO allocations for 9 assets", {
  
  skip_on_cran()

  n <- 9
  test_asset_returns <- generate_test_asset_returns(n)$returns
  if (interactive()) print(test_asset_returns)
  test_asset_correlations <- generate_test_asset_returns(n)$correlations
  if (interactive()) print(test_asset_correlations)

  optimal_allocations <- 
    calc_optimal_portfolio(
      risk_tolerance      = 0.35,
      expected_returns    = test_asset_returns$expected_return,
      standard_deviations = test_asset_returns$standard_deviation,
      correlations        = test_asset_correlations
    )$allocations  
  if (interactive()) print(optimal_allocations)

  expect_equal(
    round(optimal_allocations$total, 3),
    c(0.247, 0.221, 0.26, 0.187, 0, 0, 0, 0.084, 0)
  )
  expect_equal(sum(optimal_allocations), 1)
})

test_that("calculating optimal joint portfolio allocations", {

  skip_on_cran()
  skip_on_os("mac")
  
  n <- 9
  test_asset_returns <- generate_test_asset_returns(n)$returns
  if (interactive()) print(test_asset_returns)
  test_asset_correlations <- generate_test_asset_returns(n)$correlations
  if (interactive()) print(test_asset_correlations)
    
  effective_tax_rates <- calc_effective_tax_rate(
    portfolio                = test_asset_returns,
    tax_rate_ltcg            = 0.20,
    tax_rate_ordinary_income = 0.40
  )
  effective_tax_rates <- effective_tax_rates$aftertax$effective_tax_rate
  effective_tax_rates[7] <- 0
  expect_equal(
    round(effective_tax_rates, 4),
    c(0.2571, 0.2653, 0.2683, 0.2706, rep(0.4000, 2), 0, rep(0.4000, 2))
  )

  in_taxable_accounts      <- 250000
  in_taxadvantage_accounts <- 20500
  in_accounts <- in_taxable_accounts + in_taxadvantage_accounts
  fraction_in_taxable_accounts <- in_taxable_accounts / in_accounts
  if (interactive()) print(fraction_in_taxable_accounts)
  
  risk_tolerance <- 0.35
  
  optimal_joint_portfolio <- calc_optimal_portfolio(
    risk_tolerance      = risk_tolerance,
    expected_returns    = test_asset_returns$expected_return,
    standard_deviations = test_asset_returns$standard_deviation,
    correlations        = test_asset_correlations,
    effective_tax_rates = effective_tax_rates,
    in_taxable_accounts = fraction_in_taxable_accounts
  )$allocations 
  expect_equal(
    ignore_attr = TRUE,
    round(optimal_joint_portfolio$taxable, 3),
    c(0.271, 0.248, 0.292, 0.112, 0, 0, 0, 0, 0)
  )
  expect_equal(
    ignore_attr = TRUE,
    round(optimal_joint_portfolio$taxadvantaged, 3),
    c(0, 0, 0, 0.076, 0, 0, 0, 0, 0)
  )

  expect_equal(
    sum(optimal_joint_portfolio$taxable) + 
      sum(optimal_joint_portfolio$taxadvantaged), 
    1
  )
})

test_that("calculating optimal joint net-worth portfolio allocations", {
  
  skip_on_cran()
  skip_on_os("mac")
  
  n <- 9
  test_asset_returns <- generate_test_asset_returns(n)$returns
  if (interactive()) print(test_asset_returns)
  test_asset_correlations <- generate_test_asset_returns(n)$correlations
  if (interactive()) print(test_asset_correlations)
    
  effective_tax_rates <- calc_effective_tax_rate(
    portfolio               = test_asset_returns,
    tax_rate_ltcg            = 0.20,
    tax_rate_ordinary_income = 0.40
  )
  effective_tax_rates <- effective_tax_rates$aftertax$effective_tax_rate
  effective_tax_rates[7] <- 0
  expect_equal(
    round(effective_tax_rates, 4),
    c(0.2571, 0.2653, 0.2683, 0.2706, rep(0.4000, 2), 0, rep(0.4000, 2))
  )

  in_taxable_accounts      <- 250000
  in_taxadvantage_accounts <- 20500
  in_accounts <- in_taxable_accounts + in_taxadvantage_accounts
  fraction_in_taxable_accounts <- in_taxable_accounts / in_accounts
  fraction_in_taxable_accounts
  
  risk_tolerance <- 0.35
    
  weights <- 
    test_asset_returns |> 
    dplyr::select(name) |> 
    dplyr::mutate(
      human_capital = c(0.0936, 0.0468, 0.0468, 0, 0.3746, 0.3746, 0, 0, 0.0635),
      liabilities = c(0.1263, 0, 0, 0, 0.3368, 0.3789, 0, 0, 0.1581)
    )
  expect_equal(round(sum(weights$human_capital), 3), 1)
  expect_equal(round(sum(weights$liabilities), 3), 1)

  financial_wealth <- in_accounts
  human_capital <- 2767689
  liabilities   <- 1392064
  nondiscretionary_consumption <- 40000
  consumption <- 86000
  income <- 75000
  life_insurance_premium <- (1- 0.99999)/(1 + 0.025)*1200000

  changed_returns <- test_asset_returns$expected_return
  changed_returns[7] <- -.10

  optimal_joint_networth_portfolio = calc_optimal_portfolio(
    risk_tolerance = 0.35,
    expected_returns = changed_returns,
    standard_deviations = test_asset_returns$standard_deviation,
    correlations = test_asset_correlations,
    effective_tax_rates = effective_tax_rates,
    in_taxable_accounts = fraction_in_taxable_accounts,
    financial_wealth = financial_wealth,
    human_capital = human_capital,
    liabilities = liabilities,
    nondiscretionary_consumption = nondiscretionary_consumption,
    discretionary_consumption = consumption - nondiscretionary_consumption,
    income = income,
    life_insurance_premium = life_insurance_premium,
    human_capital_weights = weights$human_capital,
    liabilities_weights = weights$liabilities
  )$allocations  
  if (interactive())
    optimal_joint_networth_portfolio |> 
    purrr::map(function(x) {
      x |> format_percent()
    })
  
  expect_equal(
    optimal_joint_networth_portfolio |> 
      dplyr::select(taxable, taxadvantaged) |>
      purrr::map(function(x) sum(x)) |> 
      as.numeric() |> 
      sum(), 
    1
  )
  expect_equal(
    optimal_joint_networth_portfolio |> 
      purrr::map(function(x) sum(x)) |> 
      as.numeric() |> 
      round(3),
    c(0.924, 0.076, 1)
  )
  expect_equal(
    round(optimal_joint_networth_portfolio$taxable, 3),
    c(0.133, 0.021, 0.151, 0.620, rep(0, 5))
  )
  expect_equal(
    round(optimal_joint_networth_portfolio$taxadvantaged, 3),
    c(0, 0, 0.0, 0.076, rep(0, 4), 0)
  )
})
