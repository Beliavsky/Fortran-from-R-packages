test_that("Input",{
  expect_error(equilibrium_mean())
  expect_error(equilibrium_mean(risk="VAR"))
  expect_error(equilibrium_mean(dat, market_portfolio))
})
