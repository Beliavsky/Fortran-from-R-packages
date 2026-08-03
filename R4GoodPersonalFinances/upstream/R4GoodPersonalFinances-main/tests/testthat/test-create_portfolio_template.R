test_that("printing portfolio", {

  portfolio <- create_portfolio_template()
  expect_snapshot(
    portfolio
  )
})
