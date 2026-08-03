test_that("calculating gompertz retirement ruin probability", {

  expect_error(incomplete_gamma(a = 0, c = 0))
  
  expect_equal(
    calc_incomplete_gamma(a = 2, c = 3),
    0.19914,
    tolerance = 0.0001
  )
  
  expect_equal(
    calc_incomplete_gamma(a = 3, c = 2),
    1.3533,
    tolerance = 0.0001
  )
  
  expect_equal(
    calc_incomplete_gamma(a = -0.5, c = 1),
    0.1781477,
    tolerance = 0.0001
  )
  
  expect_equal(calc_incomplete_gamma(a = 1, c = 0), 1)
  expect_equal(calc_incomplete_gamma(a = 2, c = 0), 1)
  expect_equal(calc_incomplete_gamma(a = 3, c = 0), 2)

  expect_equal(
    calc_a(0.01, 65, 88, 10),
    17.88846,
    tolerance = 0.0001
  )

  expect_equal(
    calc_a(0.03, 65, 88, 10),
    14.408,
    tolerance = 0.0001
  )

  calc_retirement_ruin(
    age                   = 65,
    gompertz_mode         = 88,
    gompertz_dispersion   = 10,
    spending_rate         = 0.05,
    portfolio_return_mean = 0.02,
    portfolio_return_sd   = 0.20
  ) |> 
    expect_equal(0.3362243, tolerance = 0.00001)
  
  calc_retirement_ruin(
    age                   = 65,
    gompertz_mode         = 88,
    gompertz_dispersion   = 10,
    spending_rate         = 0.04,
    portfolio_return_mean = 0.02,
    portfolio_return_sd   = 0.20
  ) |> 
    expect_equal(0.2216641, tolerance = 0.00001)

  calc_retirement_ruin(
    age                   = 65,
    gompertz_mode         = 88,
    gompertz_dispersion   = 10,
    spending_rate         = 0.04,
    portfolio_return_mean = 0.02,
    portfolio_return_sd   = 0.15
  ) |> 
    expect_equal(0.1744054, tolerance = 0.00001)

})
  