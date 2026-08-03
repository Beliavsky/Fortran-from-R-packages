test_that("purchasing power with non-negative real interest rates", {
  
  calc_purchasing_power(x = 1, years = 1, real_interest_rate = 0) |> 
    expect_equal(1)
  
  calc_purchasing_power(x = 1, years = 50, real_interest_rate = 0.02) |> 
    round(1) |> 
    expect_equal(2.7)
})

test_that("purchasing power with negative real interest rates", {

  calc_purchasing_power(x = 10, years = 30, real_interest_rate = -0.02) |> 
    round(1) |> 
    expect_equal(5.5)
})

