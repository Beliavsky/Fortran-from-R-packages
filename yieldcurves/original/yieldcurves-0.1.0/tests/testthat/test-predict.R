test_that("yc_predict works with nelson_siegel", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  pred <- yc_predict(fit, c(3, 7, 15))
  expect_equal(nrow(pred), 3)
  expect_true(all(c("maturity", "rate") %in% names(pred)))
  expect_equal(pred$maturity, c(3, 7, 15))
})

test_that("yc_predict reproduces fitted values", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  pred <- yc_predict(fit, d$maturities)
  expect_equal(pred$rate, fit$fitted, tolerance = 1e-10)
})

test_that("yc_predict works with observed curve", {
  d <- make_simple_yields()
  curve <- yc_curve(d$maturities, d$rates)
  pred <- yc_predict(curve, c(3, 7, 15))
  expect_equal(nrow(pred), 3)
  # Should interpolate linearly
  expect_true(all(!is.na(pred$rate)))
})

test_that("yc_predict works with svensson", {
  d <- make_us_yields()
  fit <- yc_svensson(d$maturities, d$rates)
  pred <- yc_predict(fit, c(3, 7))
  expect_equal(nrow(pred), 2)
})

test_that("yc_predict works with cubic_spline", {
  d <- make_us_yields()
  fit <- yc_cubic_spline(d$maturities, d$rates)
  pred <- yc_predict(fit, c(3, 7))
  expect_equal(nrow(pred), 2)
})

test_that("yc_predict rejects non-curve", {
  expect_error(yc_predict(list(), c(1, 2)), "yc_curve")
})

test_that("yc_predict rejects invalid maturities", {
  d <- make_simple_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  expect_error(yc_predict(fit, c(-1, 2)), "positive")
})
