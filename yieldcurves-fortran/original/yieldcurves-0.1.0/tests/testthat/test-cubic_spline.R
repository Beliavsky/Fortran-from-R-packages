test_that("yc_cubic_spline fits exactly through data points", {
  d <- make_us_yields()
  fit <- yc_cubic_spline(d$maturities, d$rates)
  expect_s3_class(fit, "yc_curve")
  expect_equal(fit$method, "cubic_spline")
  # Spline should pass exactly through data points
  expect_equal(fit$residuals, rep(0, length(d$rates)), tolerance = 1e-10)
})

test_that("yc_cubic_spline predicts between points", {
  d <- make_simple_yields()
  fit <- yc_cubic_spline(d$maturities, d$rates)
  pred <- yc_predict(fit, c(3, 7, 15))
  expect_equal(nrow(pred), 3)
  # Predicted rates should be in reasonable range
  expect_true(all(pred$rate > 0.03 & pred$rate < 0.06))
})

test_that("yc_cubic_spline supports fmm method", {
  d <- make_us_yields()
  fit <- yc_cubic_spline(d$maturities, d$rates, method = "fmm")
  expect_equal(fit$params$spline_method, "fmm")
  expect_equal(fit$residuals, rep(0, length(d$rates)), tolerance = 1e-10)
})

test_that("yc_cubic_spline rejects fewer than 3 points", {
  expect_error(yc_cubic_spline(c(1, 2), c(0.04, 0.05)), "at least 3")
})
