test_that("Input",{
  expect_error(observ_powerexp())
  expect_error(observ_powerexp(covmat))
  expect_error(observ_powerexp(x,beta))
})