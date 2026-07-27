test_that("Input",{
  expect_error(.discrete_variance())
  expect_error(.discrete_variance(returns_coef))
})