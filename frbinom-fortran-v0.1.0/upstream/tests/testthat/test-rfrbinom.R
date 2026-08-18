test_that("fractional binomial distributions", {
  expect_equal(dfrbinom(4,10, .6, .6, .2), 0.079190056)
})
