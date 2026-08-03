test_that("plotting future spending", {
  
  skip_on_cran()
  skip_on_ci()

  skip_on_ci()

  older_member <- HouseholdMember$new(
    name       = "older",  
    birth_date = "1980-02-15"
  )  
  older_member$mode       <- 80
  older_member$dispersion <- 10

  younger_member <- HouseholdMember$new(
    name       = "younger",  
    birth_date = "1990-07-15"
  )
  younger_member$mode       <- 85
  younger_member$dispersion <- 9

  household <- Household$new()
  household$add_member(older_member)  
  household$add_member(younger_member)
  
  household$expected_income <- list(
    "income_1" = c(
      "members$older$age >= 44 & members$older$age < 46 ~ 100",
      "members$older$age >= 46 ~ 3000"
    )
  )
  household$expected_spending <- list(
    "spending1" = c(
      "members$older$age <= 65 ~ 3000"
    ),
    "spending2" = c(
      "TRUE ~ 16000"
    )
  )
  test_current_date <- "2020-07-15"
  portfolio <- generate_test_asset_returns(2)$returns

  scenario <- 
    simulate_scenario(
      household    = household,
      portfolio    = portfolio,
      current_date = test_current_date,
      monte_carlo_samples = 3,
      seeds = 1234
    )
  
  plot <- plot_future_spending(
    scenario = scenario,
    type     = "discretionary",
    period   = "monthly"
  ); if (interactive()) print(plot)
  vdiffr::expect_doppelganger("plot1", plot)

  plot <- plot_future_spending(
    scenario = scenario,
    type     = "discretionary"
  ); if (interactive()) print(plot)
  vdiffr::expect_doppelganger("plot4", plot)

  plot <- plot_future_spending(
    scenario = scenario,
    type     = "non-discretionary",
    period   = "monthly"
  ); if (interactive()) print(plot)
  vdiffr::expect_doppelganger("plot2", plot)

  plot <- plot_future_spending(
    scenario = scenario,
    type     = "non-discretionary"
  ); if (interactive()) print(plot)
  vdiffr::expect_doppelganger("plot5", plot)

  plot <- plot_future_spending(
    scenario = scenario,
    period   = "monthly"
  ); if (interactive()) print(plot)
  vdiffr::expect_doppelganger("plot3", plot)


  plot <- plot_future_spending(
    scenario = scenario,
    type     = "discretionary",
    period   = "monthly",
    y_limits = c(-2000, 2000)
  ); if (interactive()) print(plot)
  vdiffr::expect_doppelganger("plot6", plot)

  plot <- plot_future_spending(
    scenario = scenario,
    period   = "monthly",
    y_limits = c(-1000, 2000)
  ); if (interactive()) print(plot)
  vdiffr::expect_doppelganger("plot7", plot)
})
