test_that("yc_bond_duration returns correct structure", {
  result <- yc_bond_duration(face = 100, coupon_rate = 0.05, maturity = 2,
                              yield = 0.04, frequency = 2)
  expect_type(result, "list")
  expect_true(all(c("macaulay_duration", "modified_duration", "convexity",
                     "price") %in% names(result)))
})

test_that("yc_bond_duration: hand-computed 2Y 5% at 4% semi-annual", {
  # 2-year, 5% coupon, semi-annual, yield = 4%

  # 4 periods, coupon = 2.5, face = 100
  # CF: 2.5, 2.5, 2.5, 102.5
  # y_per = 0.02
  # PV factors: 1/1.02, 1/1.02^2, 1/1.02^3, 1/1.02^4
  # times: 0.5, 1.0, 1.5, 2.0
  pv1 <- 2.5 / 1.02
  pv2 <- 2.5 / 1.02^2
  pv3 <- 2.5 / 1.02^3
  pv4 <- 102.5 / 1.02^4

  price_expected <- pv1 + pv2 + pv3 + pv4
  mac_expected <- (0.5 * pv1 + 1.0 * pv2 + 1.5 * pv3 + 2.0 * pv4) / price_expected
  mod_expected <- mac_expected / 1.02

  # Convexity: sum(t*(t+0.5)*CF/(1+y/2)^i) / (P * (1+y/2)^2)
  conv_expected <- (0.5 * 1.0 * pv1 + 1.0 * 1.5 * pv2 +
                      1.5 * 2.0 * pv3 + 2.0 * 2.5 * pv4) / (price_expected * 1.02^2)

  result <- yc_bond_duration(face = 100, coupon_rate = 0.05, maturity = 2,
                              yield = 0.04, frequency = 2)

  expect_equal(result$price, price_expected, tolerance = 1e-6)
  expect_equal(result$macaulay_duration, mac_expected, tolerance = 1e-6)
  expect_equal(result$modified_duration, mod_expected, tolerance = 1e-6)
  expect_equal(result$convexity, conv_expected, tolerance = 1e-6)
})

test_that("yc_bond_duration: zero coupon is just maturity", {
  result <- yc_bond_duration(face = 100, coupon_rate = 0, maturity = 5,
                              yield = 0.04, frequency = 2)
  expect_equal(result$macaulay_duration, 5, tolerance = 1e-6)
})

test_that("yc_bond_duration: annual compounding", {
  result <- yc_bond_duration(face = 100, coupon_rate = 0.06, maturity = 3,
                              yield = 0.05, frequency = 1,
                              compounding = "annual")
  expect_true(result$macaulay_duration > 0)
  expect_true(result$modified_duration < result$macaulay_duration)
  expect_true(result$price > 100) # coupon > yield means premium bond
})

test_that("yc_bond_duration: continuous compounding", {
  result <- yc_bond_duration(face = 100, coupon_rate = 0.05, maturity = 2,
                              yield = 0.04, frequency = 2,
                              compounding = "continuous")
  expect_true(result$macaulay_duration > 0)
  # For continuous, modified = macaulay
  expect_equal(result$macaulay_duration, result$modified_duration, tolerance = 1e-10)
})

test_that("yc_bond_duration: higher coupon means lower duration", {
  d1 <- yc_bond_duration(face = 100, coupon_rate = 0.02, maturity = 10,
                          yield = 0.04, frequency = 2)
  d2 <- yc_bond_duration(face = 100, coupon_rate = 0.08, maturity = 10,
                          yield = 0.04, frequency = 2)
  expect_true(d1$macaulay_duration > d2$macaulay_duration)
})

test_that("yc_bond_duration rejects bad inputs", {
  expect_error(yc_bond_duration(face = -100, coupon_rate = 0.05, maturity = 2,
                                 yield = 0.04), "positive")
  expect_error(yc_bond_duration(face = 100, coupon_rate = "a", maturity = 2,
                                 yield = 0.04), "numeric")
  expect_error(yc_bond_duration(face = 100, coupon_rate = 0.05, maturity = 2,
                                 yield = 0.04, frequency = 4), "frequency")
})
