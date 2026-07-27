test_that("yc_carry returns correct structure", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  result <- yc_carry(fit)
  expect_true(is.data.frame(result))
  expect_true(all(c("maturity", "carry", "rolldown", "total") %in% names(result)))
})

test_that("yc_carry total equals carry plus rolldown", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  result <- yc_carry(fit)
  expect_equal(result$total, result$carry + result$rolldown, tolerance = 1e-10)
})

test_that("yc_carry at specific maturities", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  result <- yc_carry(fit, maturities = c(2, 5, 10))
  expect_equal(nrow(result), 3)
  expect_equal(result$maturity, c(2, 5, 10))
})

test_that("yc_carry custom horizon", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  result_1m <- yc_carry(fit, maturities = c(2, 10), horizon = 1 / 12)
  result_3m <- yc_carry(fit, maturities = c(2, 10), horizon = 3 / 12)
  # 3-month carry should be roughly 3x one-month carry
  expect_true(all(abs(result_3m$carry) > abs(result_1m$carry) * 2))
})

test_that("yc_carry rejects maturities <= horizon", {
  d <- make_simple_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  expect_error(yc_carry(fit, maturities = c(0.01, 2), horizon = 1 / 12), "greater than")
})

test_that("yc_carry on flat curve gives zero rolldown", {
  d <- make_flat_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  result <- yc_carry(fit, maturities = c(2, 5, 10))
  # On a flat curve, roll-down should be approximately 0
  expect_equal(result$rolldown, rep(0, 3), tolerance = 0.001)
})

test_that("yc_carry rejects non-curve", {
  expect_error(yc_carry(list()), "yc_curve")
})
