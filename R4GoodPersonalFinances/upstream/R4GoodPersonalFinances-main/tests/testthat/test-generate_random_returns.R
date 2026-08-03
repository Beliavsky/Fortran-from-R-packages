test_that("generating random returns", {

  portfolio <- generate_test_asset_returns(2)$returns
    
  returns <- 
    generate_random_returns(
      portfolio = portfolio, 
      n = 100000,
      seed = 1234
    )
  
  expect_snapshot(returns)
  expect_equal(NROW(returns), 100000)

  expect_equal(
    mean(returns[[1]]),
    portfolio$expected_return[[1]],
    tolerance = 0.01
  )
  expect_equal(
    sd(returns[[1]]),
    portfolio$standard_deviation[[1]],
    tolerance = 0.01
  )

  expect_equal(
    mean(returns[[2]]),
    portfolio$expected_return[[2]],
    tolerance = 0.01
  )
  expect_equal(
    sd(returns[[2]]),
    portfolio$standard_deviation[[2]],
    tolerance = 0.01
  )
})
