test_that("calculating present value", {

  test_cashflow <- rep(100000, 10)
  
  pv <- calc_present_value(
    cashflow      = test_cashflow,
    discount_rate = 0.02
  )

  expect_equal(
    pv[1],
    916223.671
  )
})

test_that("calculating present value", {

  test_cashflow <- c(
    rep(100000, 10),
    rep(0, 10),
    rep(10000, 20),
    rep(0, 10)
  )
  if (interactive()) format_currency(test_cashflow) 
  
  pv <- calc_present_value(
    cashflow      = test_cashflow,
    discount_rate = 0.02
  )

  expect_equal(
    pv[1],
    1028464.938
  )
  expect_snapshot(pv)
})
