test_that("plotting household survival", {
  
  skip_on_cran()
  skip_on_ci()

  test_current_date <- "2020-01-01"
  hm1 <- 
    HouseholdMember$new(
      name       = "member1",
      birth_date = "1955-01-01",
      mode       = 88,
      dispersion = 10.65
    )
  hm2 <- 
    HouseholdMember$new(
      name       = "member2",
      birth_date = "1965-01-01",
      mode       = 91,
      dispersion = 8.88
    )
  hm3 <- 
    HouseholdMember$new(
      name       = "member3",
      birth_date = "1975-01-01",
      mode       = 88,
      dispersion = 7.77
    )
  household <- Household$new()
  household$add_member(hm1)
  household$add_member(hm2)
  household$add_member(hm3) 
  
  plot <- plot_survival(
    household    = household, 
    current_date = test_current_date
  )
  if (interactive()) print(plot)
  vdiffr::expect_doppelganger("plot1", plot)
})
