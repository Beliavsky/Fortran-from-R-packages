test_that("calculate gompertz survival probability", {

  # Test based on: Milevsky, M. A., Robinson, C. (2000). Self-annuitization and ruin in retirement. North American Actuarial Journal, 4(4), p. 114.
  
  calc_gompertz_survival_probability(
    current_age = 65,
    target_age  = 85,
    mode        = 80,
    dispersion  = 10
  ) |> expect_equal(0.2404, tolerance = 0.001)

  calc_gompertz_survival_probability(
    current_age = 75,
    target_age  = 85,
    mode        = 80,
    dispersion  = 10
  ) |> expect_equal(0.3527, tolerance = 0.001)  

  calc_gompertz_survival_probability(
    current_age = 65,
    target_age  = 85,
    mode        = 81.95,
    dispersion  = 10.6
  ) |> expect_equal(0.3226, tolerance = 0.001)  

  calc_gompertz_survival_probability(
    current_age = 65,
    target_age  = 85,
    mode        = 87.8,
    dispersion  = 9.5
  ) |> expect_equal(0.5199, tolerance = 0.001)  

})

test_that("calibrating gompertz model for males on test data", {

  mortality_rates <- 
    test_mortality_rates |> 
    dplyr::filter(sex == "male")

  expect_equal(tail(mortality_rates$mortality_rate, 1), 1)
  
  params <- 
    mortality_rates |> 
      calc_gompertz_parameters(current_age = 65)
    
  expect_equal(head(params$data$survival_rate, 1), 1)
  expect_equal(tail(params$data$survival_rate, 1), 0)
  expect_equal(params$data$probability_of_death |> sum(), 1)
  
  pod_males <- 
    function() {
      plot(
        x = params$data$age,
        y = params$data$probability_of_death
      )
    }
  if (interactive()) print(pod_males())
  vdiffr::expect_doppelganger("pod_males", pod_males)
  
  expect_equal(params$mode, 86)
  expect_equal(params$dispersion, 10.48, tolerance = 0.01)

  gc_males <- function() plot_gompertz_calibration(params = params)
  if (interactive()) print(gc_males())
  vdiffr::expect_doppelganger("gc_males", gc_males)
})

test_that("calibrating gompertz model for females on test data", {
  
  mortality_rates <- 
    test_mortality_rates |> 
    dplyr::filter(sex == "female")
  
  expect_equal(tail(mortality_rates$mortality_rate, 1), 1)
  
  params <- 
    mortality_rates |> 
      calc_gompertz_parameters(current_age = 65)
    
  expect_equal(head(params$data$survival_rate, 1), 1)
  expect_equal(tail(params$data$survival_rate, 1), 0)
  expect_equal(params$data$probability_of_death |> sum(), 1)

  pod_females <- 
    function() {
      plot(
        x = params$data$age,
        y = params$data$probability_of_death
      )
    }
  if (interactive()) print(pod_females())
  vdiffr::expect_doppelganger("pod_females", pod_females)

  expect_equal(params$mode, 90)
  expect_equal(params$dispersion, 8.63, tolerance = 0.01)
  
  gc_females <- function() plot_gompertz_calibration(params = params)
  if (interactive()) print(gc_females())
  vdiffr::expect_doppelganger("gc_females", gc_females)
})

test_that("calibrating gompertz model on HMD data", {
  
  mortality_rates <- 
    life_tables |> 
    dplyr::filter(
      country == "USA" & 
      sex     == "female"
    ) |>
    dplyr::filter(year == max(year)) 
    
  params <- 
    mortality_rates |> 
      calc_gompertz_parameters(
    estimate_max_age = TRUE,
    current_age      = 0
  ) 
  
  hmd <- function() plot_gompertz_calibration(params = params)
  if (interactive()) print(hmd())
  vdiffr::expect_doppelganger("hmd", hmd)

})

test_that("calibrating joint gompertz model", {

  params <- 
    calc_gompertz_joint_parameters(
      p1 = list(
        age        = 65,
        mode       = 88,
        dispersion = 10.65
      ),
      p2 = list(
        age        = 65,
        mode       = 91,
        dispersion = 8.88
      ),
      max_age = 110
    )
  if (interactive()) print(params)

  expect_equal(
    params$mode, 
    94.764246
  )
  expect_equal(
    params$dispersion, 
    6.18133077
  )
  
  js <- function() plot_joint_survival(params = params)
  if (interactive()) print(js())
  vdiffr::expect_doppelganger("js", js)

  jsg <- function() plot_joint_survival(
    params = params, 
    include_gompertz = TRUE
  )
  if (interactive()) print(jsg())
  vdiffr::expect_doppelganger("jsg", jsg)
})

test_that("survival probability is not NaN", {

  calc_gompertz_survival_probability(
    current_age = 100,
    target_age  = 100,
    mode        = 80,
    dispersion  = 10,
    max_age     = 100
  ) |> expect_equal(0)
})
