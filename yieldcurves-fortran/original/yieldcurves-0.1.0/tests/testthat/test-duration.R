test_that("yc_duration returns correct structure", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  result <- yc_duration(fit)
  expect_true(is.data.frame(result))
  expect_true(all(c("maturity", "macaulay_duration", "modified_duration",
                     "convexity") %in% names(result)))
  expect_equal(nrow(result), length(d$maturities))
})

test_that("yc_duration continuous: duration = maturity for zero-coupon", {
  curve <- yc_curve(c(1, 5, 10), c(0.05, 0.04, 0.04))
  result <- yc_duration(curve, compounding = "continuous")
  expect_equal(result$macaulay_duration, c(1, 5, 10))
  expect_equal(result$modified_duration, c(1, 5, 10))
})

test_that("yc_duration continuous: convexity = maturity^2", {
  curve <- yc_curve(c(1, 5, 10), c(0.05, 0.04, 0.04))
  result <- yc_duration(curve, compounding = "continuous")
  expect_equal(result$convexity, c(1, 25, 100))
})

test_that("yc_duration annual compounding", {
  curve <- yc_curve(c(1, 10), c(0.05, 0.05))
  result <- yc_duration(curve, compounding = "annual")
  # Modified duration = maturity / (1 + r)
  expect_equal(result$modified_duration[1], 1 / 1.05, tolerance = 1e-10)
  expect_equal(result$modified_duration[2], 10 / 1.05, tolerance = 1e-10)
})

test_that("yc_duration semi-annual compounding", {
  curve <- yc_curve(c(1, 10), c(0.06, 0.06))
  result <- yc_duration(curve, compounding = "semi_annual")
  # Modified duration = maturity / (1 + r/2)
  expect_equal(result$modified_duration[1], 1 / 1.03, tolerance = 1e-10)
  expect_equal(result$modified_duration[2], 10 / 1.03, tolerance = 1e-10)
})

test_that("yc_duration at custom maturities", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  result <- yc_duration(fit, maturities = c(2, 10, 30))
  expect_equal(nrow(result), 3)
  expect_equal(result$maturity, c(2, 10, 30))
})

test_that("yc_duration modified < macaulay for annual compounding", {
  curve <- yc_curve(c(1, 5, 10), c(0.05, 0.04, 0.04))
  result <- yc_duration(curve, compounding = "annual")
  expect_true(all(result$modified_duration < result$macaulay_duration))
})

test_that("yc_duration rejects non-curve", {
  expect_error(yc_duration(list()), "yc_curve")
})
