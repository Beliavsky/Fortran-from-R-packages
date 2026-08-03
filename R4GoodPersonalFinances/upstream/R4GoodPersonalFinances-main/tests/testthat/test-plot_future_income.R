test_that("plotting future income", {
  
  skip_on_cran()
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
      "members$older$age <= 65 ~ 1000 * 12",
      "members$older$age >= 80 ~ 4000"
    ),
    "income_2" = c(
      "members$older$age <= 65 ~ 3000 * 12",
      "members$older$age >= 70 ~ 2000 * 12"
    )
  )
  household$expected_spending <- list(
    "spending1" = c(
      "members$older$age <= 65 ~ 3000 * 12"
    ),
    "spending2" = c(
      "TRUE ~ 2000 * 12"
    )
  )
  test_current_date <- "2020-07-15"
  portfolio <- generate_test_asset_returns(2)$returns

  scenario <- 
    simulate_scenario(
      household    = household,
      portfolio    = portfolio,
      current_date = test_current_date
    )
  
  plot <- plot_future_income(
    scenario = scenario,
    period   = "monthly"
  ); if (interactive()) print(plot)
  vdiffr::expect_doppelganger("plot1", plot)
  
  plot <- plot_future_income(
    scenario = scenario
  ); if (interactive()) print(plot)
  vdiffr::expect_doppelganger("plot2", plot)

})