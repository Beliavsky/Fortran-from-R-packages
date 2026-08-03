test_that("calculating life expectancy", {
  
  life_expectancy <- 
    calc_life_expectancy(
      current_age = 25,
      mode        = 91,
      dispersion  = 8.88,
      max_age     = 115
    )
  if (interactive()) print(life_expectancy)
  
  expect_equal(round(life_expectancy), 86)
})

test_that("calculating mode for new life expectancy", {

  mode <- 
    calc_gompertz_mode(
      life_expectancy = 86,
      current_age     = 25,
      dispersion      = 8.88,
      max_age         = 115
    )
  if (interactive()) print(mode)
  expect_equal(round(mode), 91)

  mode <- 
    calc_gompertz_mode(
      life_expectancy = 94,
      current_age     = 25,
      dispersion      = 8.88,
      max_age         = 115
    )
  if (interactive()) print(mode)
  expect_equal(round(mode), 99)

  life_expectancy <- 
    calc_life_expectancy(
      current_age = 25,
      mode        = 99,
      dispersion  = 8.88
    )
  expect_equal(round(life_expectancy), 94)
})
