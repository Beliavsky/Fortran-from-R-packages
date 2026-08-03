test_that("plotting life expectancy", {
  
  skip_on_cran()
  skip_on_ci()

  params <- list(
    default = list(mode = 91,       dispersion = 8.88),
    isabela = list(mode = 99.16564, dispersion = 8.88)
  )
  
  plot <- plot_life_expectancy(
    household = params
  ); if (interactive()) print(plot)
  vdiffr::expect_doppelganger("plot_lep", plot)
})

test_that("plotting life expectancy of a household", {
  
  skip_on_cran()
  skip_on_ci()

  hm1 <- 
    HouseholdMember$new(
      name       = "isabela",  
      birth_date = Sys.Date() - lubridate::years(25),
      mode       = 99.10791, 
      dispersion = 8.88
    )  

  household <- Household$new()
  household$add_member(hm1)  
  
  plot <- plot_life_expectancy(
    household = household
  ); if (interactive()) print(plot)
  vdiffr::expect_doppelganger("plot_leh", plot)

  hm1$calc_life_expectancy(current_date = Sys.Date()) |> 
    expect_equal(94, tolerance = 1e-5)
})
