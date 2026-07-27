test_that("yc_discount computes continuous discount factors", {
  curve <- yc_curve(c(1, 2, 5), c(0.05, 0.05, 0.05))
  df <- yc_discount(curve)
  expect_equal(df$discount_factor, exp(-0.05 * c(1, 2, 5)), tolerance = 1e-10)
})

test_that("yc_discount computes annual discount factors", {
  curve <- yc_curve(c(1, 2, 5), c(0.05, 0.05, 0.05))
  df <- yc_discount(curve, compounding = "annual")
  expect_equal(df$discount_factor, (1.05)^(-c(1, 2, 5)), tolerance = 1e-10)
})

test_that("yc_discount computes semi-annual discount factors", {
  curve <- yc_curve(c(1, 2), c(0.05, 0.05))
  df <- yc_discount(curve, compounding = "semi_annual")
  expect_equal(df$discount_factor, (1 + 0.05 / 2)^(-2 * c(1, 2)), tolerance = 1e-10)
})

test_that("yc_discount at custom maturities", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  df <- yc_discount(fit, maturities = c(1, 5, 10))
  expect_equal(nrow(df), 3)
  expect_true(all(df$discount_factor > 0 & df$discount_factor < 1))
})

test_that("yc_discount at maturity 0 gives 1 (approximately)", {
  # Very short maturity should have DF close to 1
  curve <- yc_curve(c(0.01, 1, 2), c(0.05, 0.05, 0.05))
  df <- yc_discount(curve, maturities = 0.001)
  expect_equal(df$discount_factor, 1, tolerance = 0.001)
})

test_that("yc_discount rejects non-curve", {
  expect_error(yc_discount(list()), "yc_curve")
})
