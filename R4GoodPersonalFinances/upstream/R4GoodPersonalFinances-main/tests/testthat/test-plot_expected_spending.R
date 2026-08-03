test_that("plotting expected spending", {
  
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
  
  household$expected_income <- list(
    "income_older" = c(
      "members$older$age >= 44 & members$older$age < 55 ~ 1000",
      "members$older$age >= 55 ~ 3000"
    )
  )
  household$expected_spending <- list(
    "spending1" = c(
      "members$older$age >= 65 ~ 3000",
      "TRUE ~ 6000"
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
  
  es1 <- function() plot_expected_spending(scenario)
  if (interactive()) print(es1())
  vdiffr::expect_doppelganger("es1", es1)

  es2 <- function() plot_expected_spending(
    scenario, 
    discretionary_spending_position = "bottom"
  )
  if (interactive()) print(es2())
  vdiffr::expect_doppelganger("es2", es2)

  es3 <- function() plot_expected_spending(
    scenario, 
    period = "monthly"
  )
  if (interactive()) print(es3())
  vdiffr::expect_doppelganger("es3", es3)
})
