test_that("calculating optimal asset allocation", {
  
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
      "members$older$age < 60 ~ 40000",
      "members$older$age >= 60 ~ 20000"
    )
  )
  household$expected_spending <- list(
    "spending1" = c(
      "TRUE ~ 10000"
    )
  )
  current_date <- "2020-07-15"

  portfolio <- generate_test_asset_returns(2)$returns

  optimal_allocation <- 
    calc_optimal_asset_allocation(
      household    = household,
      portfolio    = portfolio,
      current_date = current_date
    )
  
  if (interactive()) print(optimal_allocation, width = Inf)
  
  expect_true(tibble::is_tibble(optimal_allocation$allocations))
  expect_true(tibble::is_tibble(optimal_allocation$allocations$optimal))
  expect_true(tibble::is_tibble(optimal_allocation$allocations$current))
  
  expect_true("total" %in% names(optimal_allocation$allocations$optimal))
  expect_true("total" %in% names(optimal_allocation$allocations$current))
  expect_true("taxable" %in% names(optimal_allocation$allocations$optimal))
  expect_true("taxable" %in% names(optimal_allocation$allocations$current))
  expect_true(
    "taxadvantaged" %in% names(optimal_allocation$allocations$current)
  )
  expect_true(
    "taxadvantaged" %in% names(optimal_allocation$allocations$optimal)
  )
  expect_equal(
    sum(optimal_allocation$allocations$optimal$total),
    1
  )
  expect_equal(
    sum(optimal_allocation$allocations$current$total),
    1
  )
})
