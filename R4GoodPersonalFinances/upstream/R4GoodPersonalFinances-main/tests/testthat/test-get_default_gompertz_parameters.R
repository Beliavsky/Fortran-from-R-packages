test_that("getting default gompertz parameters", {
    
  gompertz_parameters <- 
    get_default_gompertz_parameters(
      age     = 22, 
      country = "USA", 
      sex     = "female"
    )
  
  expect_true(is.list(gompertz_parameters))
  expect_true("mode"        %in% names(gompertz_parameters))
  expect_true("dispersion"  %in% names(gompertz_parameters))
  expect_true("current_age" %in% names(gompertz_parameters))
  expect_true("max_age"     %in% names(gompertz_parameters))

  gompertz_parameters_22_usa_f <- gompertz_parameters
  
  gompertz_parameters <- 
    get_default_gompertz_parameters(
      age     = 65, 
      country = "Poland", 
      sex     = "male"
    )
  
  expect_true(
    gompertz_parameters$mode < gompertz_parameters_22_usa_f$mode
  )
  expect_true(
    gompertz_parameters$dispersion > gompertz_parameters_22_usa_f$dispersion
  )
})
