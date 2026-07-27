test_that("Input",{
  expect_error(BL_post_distr())
  expect_error(BL_post_distr(dat,market_portfolio))
})
