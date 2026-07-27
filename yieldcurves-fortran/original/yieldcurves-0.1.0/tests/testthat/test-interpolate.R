test_that("yc_interpolate linear works", {
  curve <- yc_curve(c(1, 2, 10), c(0.04, 0.05, 0.06))
  result <- yc_interpolate(curve, c(1.5, 5))
  expect_equal(nrow(result), 2)
  # Linear interpolation: 1.5 should be midpoint of 0.04 and 0.05
  expect_equal(result$rate[1], 0.045, tolerance = 1e-10)
})

test_that("yc_interpolate log_linear works", {
  curve <- yc_curve(c(1, 10), c(0.04, 0.06))
  result <- yc_interpolate(curve, 5, method = "log_linear")
  expect_equal(nrow(result), 1)
  expect_true(result$rate > 0.04 & result$rate < 0.06)
})

test_that("yc_interpolate cubic works", {
  d <- make_us_yields()
  curve <- yc_curve(d$maturities, d$rates)
  result <- yc_interpolate(curve, c(3, 7, 15), method = "cubic")
  expect_equal(nrow(result), 3)
})

test_that("yc_interpolate uses predict for fitted curves", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  interp <- yc_interpolate(fit, c(3, 7))
  pred <- yc_predict(fit, c(3, 7))
  expect_equal(interp$rate, pred$rate)
})

test_that("yc_interpolate rejects bad inputs", {
  expect_error(yc_interpolate(list(), c(1, 2)), "yc_curve")
  d <- make_simple_yields()
  curve <- yc_curve(d$maturities, d$rates)
  expect_error(yc_interpolate(curve, c(-1)), "positive")
})
