test_that("plotting scenarios metrics without Monte Carlo samples", {

  skip_on_cran()
  skip_on_ci()

  older_member <- HouseholdMember$new(
    name       = "older",  
    birth_date = "1980-02-15"
  )  
  older_member$mode       <- 80
  older_member$dispersion <- 10
  older_member$set_event("retirement", 65)

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
    "income_older" = c(
      "!members$older$events$retirement$on ~ 10000 * 12",
      "members$older$events$retirement$on ~ 1000 * 12"
    )
  )
  household$expected_spending <- list(
    "spending1" = c(
      "TRUE ~ 6000 * 12"
    )
  )
  
  portfolio <- generate_test_asset_returns(2)$returns

  start_ages <- c(55, 65, 75, 88)

  scenarios_parameters <- 
    tibble::tibble(
      member    = "older",
      event      = "retirement",
      start_age = start_ages,
      years     = Inf,
      end_age   = Inf
    ) |> 
    dplyr::mutate(
      scenario_id = start_age
    ) |> 
    tidyr::nest(events = -scenario_id)
  
  test_current_date <- "2020-07-15"

  scenarios <- 
    simulate_scenarios(
      scenarios_parameters = scenarios_parameters,
      household            = household,
      portfolio            = portfolio,
      current_date         = test_current_date
    )
  
  # set.seed(1234)
  withr::local_seed(1234)
  
    plot <- plot_scenarios(
      scenarios = scenarios
    ); if (interactive()) print(plot)
    vdiffr::expect_doppelganger("plot1", plot)
  
    plot <- plot_scenarios(
      scenarios = scenarios,
      period    = "monthly"
    ); if (interactive()) print(plot)
    vdiffr::expect_doppelganger("plot2", plot)
})

test_that("plotting scenarios metrics with Monte Carlo samples", {

  skip_on_cran()
  skip_on_ci()

  older_member <- HouseholdMember$new(
    name       = "older",  
    birth_date = "1980-02-15"
  )  
  older_member$mode       <- 80
  older_member$dispersion <- 10
  older_member$set_event("retirement", 65)

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
    "income_older" = c(
      "!members$older$events$retirement$on ~ 10000 * 12",
      "members$older$events$retirement$on ~ 1000 * 12"
    )
  )
  household$expected_spending <- list(
    "spending1" = c(
      "TRUE ~ 6000 * 12"
    )
  )
  
  portfolio <- generate_test_asset_returns(2)$returns

  start_ages <- c(45, 50, 55, 65, 75, 88)

  scenarios_parameters <- 
    tibble::tibble(
      member    = "older",
      event      = "retirement",
      start_age = start_ages,
      years     = Inf,
      end_age   = Inf
    ) |> 
    dplyr::mutate(
      scenario_id = start_age
    ) |> 
    tidyr::nest(events = -scenario_id)
  
  test_current_date <- "2020-07-15"

  scenarios <- 
    simulate_scenarios(
      scenarios_parameters = scenarios_parameters,
      household            = household,
      portfolio            = portfolio,
      current_date         = test_current_date,
      monte_carlo_samples  = 1,
      seeds = 1234
    )
  
  withr::local_seed(1234)

  plot <- plot_scenarios(
    scenarios = scenarios
  ); if (interactive()) print(plot)
  vdiffr::expect_doppelganger("plot_smmc", plot)

  plot <- plot_scenarios(
    scenarios = scenarios,
    period    = "monthly"
  ); if (interactive()) print(plot)
  vdiffr::expect_doppelganger("plot_smmcm", plot)
})
