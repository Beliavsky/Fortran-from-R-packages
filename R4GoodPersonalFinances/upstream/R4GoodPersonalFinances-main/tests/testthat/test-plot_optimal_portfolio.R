test_that("plotting optimal portfolio", {
  
  skip_on_cran()
  skip_on_ci()
  
  older_member <- HouseholdMember$new(
    name       = "older",  
    birth_date = "1980-02-15"
  )  
  older_member$mode       <- 80
  older_member$dispersion <- 10

  household <- Household$new()
  household$add_member(older_member)  
  household$risk_tolerance <- 0.1
  
  household$expected_income <- list(
    "income_older" = c(
    )
  )
  household$expected_spending <- list(
    "spending1" = c(
    )
  )
  current_date <- "2020-07-15"

  portfolio <- generate_test_asset_returns(2)$returns
  portfolio$accounts$taxable <- 
    portfolio$accounts$taxable + 0.66
  if (interactive()) print(portfolio, width = Inf)

  optimal_portfolio <- 
    calc_optimal_asset_allocation(
      household    = household,
      portfolio    = portfolio,
      current_date = current_date
    )
  
  plot <- plot_optimal_portfolio(
    portfolio = optimal_portfolio
  ); if (interactive()) print(plot)
  vdiffr::expect_doppelganger("plot1", plot)
})
