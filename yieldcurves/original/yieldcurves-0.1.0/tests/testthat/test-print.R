test_that("print.yc_curve returns invisibly", {
  d <- make_simple_yields()
  curve <- yc_curve(d$maturities, d$rates)
  result <- withVisible(print(curve))
  expect_false(result$visible)
})

test_that("print.yc_curve runs without error for NS", {
  d <- make_simple_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  expect_no_error(print(fit))
})

test_that("print.yc_curve runs without error with date", {
  curve <- yc_curve(c(1, 2), c(0.04, 0.05), date = as.Date("2024-01-15"))
  expect_no_error(print(curve))
})

test_that("print.yc_curve runs without error for fitted", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  expect_no_error(print(fit))
})

test_that("summary.yc_curve returns invisibly", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  result <- withVisible(summary(fit))
  expect_false(result$visible)
})

test_that("plot.yc_curve runs without error", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  expect_no_error(plot(fit))
})

test_that("plot.yc_curve works for observed curve", {
  d <- make_simple_yields()
  curve <- yc_curve(d$maturities, d$rates)
  expect_no_error(plot(curve))
})
