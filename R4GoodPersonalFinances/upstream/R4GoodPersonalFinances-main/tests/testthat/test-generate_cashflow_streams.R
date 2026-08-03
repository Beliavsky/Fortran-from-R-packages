test_that("generating cashflow streams", {

  h <- Household$new()
  h$add_member(
    HouseholdMember$new(
      name       = "older",  
      birth_date = "1980-02-15",
      mode       = 80,
      dispersion = 10
    )  
  )  
  h$add_member(
    HouseholdMember$new(
      name       = "younger",  
      birth_date = "1990-07-15",
      mode       = 85,
      dispersion = 9
    )
  )  

  test_current_date <- "2020-07-15"

  timeline <- 
    generate_household_timeline(
      household    = h, 
      current_date = test_current_date
    ) 
  
  test_triggers <- list(
    "income1" = c(
      "members$older$age >= 44 & members$older$age < 46 ~ 100",
      "members$older$age >= 46 ~ 300"
    ),
    "income2" = c(
      "members$younger$age >= 34 & members$younger$age < 36 ~ 44",
      "members$younger$age >= 36 ~ 55"
    )
  )
  
  cashflows <- generate_cashflow_streams(
    timeline = timeline,
    triggers = test_triggers
  ) 
  if (interactive()) 
    timeline |> 
      dplyr::bind_cols(cashflows) |>
      print(width = Inf, n = Inf)
  
  expect_equal(
    sum(cashflows$income1),
    16700
  )
  expect_equal(
    sum(cashflows$income2),
    3663
  )
})
