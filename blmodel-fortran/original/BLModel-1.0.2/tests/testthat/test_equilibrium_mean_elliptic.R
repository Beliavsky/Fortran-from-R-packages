test_that("Input",{
  expect_error(.equilibrium_mean_elliptic())
  expect_error(.equilibrium_mean_elliptic(MCov, MarketPortfolio))
})