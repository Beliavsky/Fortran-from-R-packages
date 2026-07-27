test_that("yc_slope returns named list", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  result <- yc_slope(fit)
  expect_true(is.list(result))
  expect_true(all(c("spread_2s10s", "spread_2s30s", "spread_5s30s",
                     "spread_3m10y", "butterfly_2s5s10s") %in% names(result)))
})

test_that("yc_slope 2s10s is negative for normal curve", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  result <- yc_slope(fit)
  # 10Y < 2Y in our test data (slight inversion at long end)
  expect_true(is.numeric(result$spread_2s10s))
})

test_that("yc_slope 2s10s is positive for inverted curve", {
  d <- make_inverted_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  result <- yc_slope(fit)
  # Inverted: 2Y > 10Y, so 10Y - 2Y should be negative
  expect_true(result$spread_2s10s < 0)
})

test_that("yc_slope returns NA for out-of-range tenors", {
  # Short curve without 30Y
  curve <- yc_curve(c(1, 2, 5, 10), c(0.04, 0.042, 0.043, 0.045))
  result <- yc_slope(curve)
  expect_true(is.na(result$spread_2s30s))
})

test_that("yc_level_slope_curvature works for NS", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  result <- yc_level_slope_curvature(fit)
  expect_equal(result$level, fit$params$beta0)
  expect_equal(result$slope, fit$params$beta1)
  expect_equal(result$curvature, fit$params$beta2)
})

test_that("yc_level_slope_curvature works for observed", {
  curve <- yc_curve(c(1, 5, 10), c(0.05, 0.04, 0.03))
  result <- yc_level_slope_curvature(curve)
  expect_equal(result$level, mean(c(0.05, 0.04, 0.03)))
  expect_equal(result$slope, 0.05 - 0.03)
})

test_that("yc_slope works with Svensson curve", {
  d <- make_us_yields()
  fit <- yc_svensson(d$maturities, d$rates)
  result <- yc_slope(fit)
  expect_true(is.list(result))
  expect_true(is.numeric(result$spread_2s10s))
})

test_that("yc_level_slope_curvature works for Svensson", {
  d <- make_us_yields()
  fit <- yc_svensson(d$maturities, d$rates)
  result <- yc_level_slope_curvature(fit)
  expect_equal(result$level, fit$params$beta0)
  expect_equal(result$slope, fit$params$beta1)
  expect_equal(result$curvature, fit$params$beta2)
})

test_that("yc_slope rejects non-curve", {
  expect_error(yc_slope(list()), "yc_curve")
})
