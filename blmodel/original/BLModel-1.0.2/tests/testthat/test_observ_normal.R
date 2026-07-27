test_that("Input",{
  expect_error(observ_normal())
  expect_error(observ_normal(p))
  expect_error(observ_normal(p,q))
  expect_error(observ_normal(x,p,covmat))
  })