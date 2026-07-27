test_that("Input",{
  expect_error(observ_ts())
  expect_error(observ_ts(covmat,df))
  expect_error(observ_ts(df))
})