test_that("yc_zspread: zero spread when price equals benchmark PV", {
  curve <- yc_curve(c(1, 2, 5, 10), c(0.04, 0.042, 0.044, 0.045))

  # Price the bond at the benchmark curve (no spread)
  maturity <- 5
  freq <- 2
  step <- 1 / freq
  cf_times <- seq(step, maturity, by = step)
  coupon <- 100 * 0.04 / freq
  cfs <- rep(coupon, length(cf_times))
  cfs[length(cfs)] <- cfs[length(cfs)] + 100

  z <- approx(c(1, 2, 5, 10), c(0.04, 0.042, 0.044, 0.045),
              xout = cf_times, rule = 2)$y
  benchmark_price <- sum(cfs * (1 + z)^(-cf_times))

  result <- yc_zspread(price = benchmark_price, coupon_rate = 0.04,
                        maturity = 5, curve = curve, frequency = 2)
  expect_equal(result$zspread, 0, tolerance = 1e-6)
})

test_that("yc_zspread: positive spread when price is below benchmark", {
  curve <- yc_curve(c(1, 2, 5, 10), c(0.04, 0.042, 0.044, 0.045))

  # Price the bond at the benchmark curve
  maturity <- 5
  freq <- 2
  step <- 1 / freq
  cf_times <- seq(step, maturity, by = step)
  coupon <- 100 * 0.04 / freq
  cfs <- rep(coupon, length(cf_times))
  cfs[length(cfs)] <- cfs[length(cfs)] + 100

  z <- approx(c(1, 2, 5, 10), c(0.04, 0.042, 0.044, 0.045),
              xout = cf_times, rule = 2)$y
  benchmark_price <- sum(cfs * (1 + z)^(-cf_times))

  # Lower price means positive spread
  result <- yc_zspread(price = benchmark_price - 2, coupon_rate = 0.04,
                        maturity = 5, curve = curve, frequency = 2)
  expect_true(result$zspread > 0)
})

test_that("yc_zspread: negative spread when price is above benchmark", {
  curve <- yc_curve(c(1, 2, 5, 10), c(0.04, 0.042, 0.044, 0.045))

  maturity <- 5
  freq <- 2
  step <- 1 / freq
  cf_times <- seq(step, maturity, by = step)
  coupon <- 100 * 0.04 / freq
  cfs <- rep(coupon, length(cf_times))
  cfs[length(cfs)] <- cfs[length(cfs)] + 100

  z <- approx(c(1, 2, 5, 10), c(0.04, 0.042, 0.044, 0.045),
              xout = cf_times, rule = 2)$y
  benchmark_price <- sum(cfs * (1 + z)^(-cf_times))

  result <- yc_zspread(price = benchmark_price + 1, coupon_rate = 0.04,
                        maturity = 5, curve = curve, frequency = 2)
  expect_true(result$zspread < 0)
})

test_that("yc_zspread: model_price matches input price", {
  curve <- yc_curve(c(1, 2, 5, 10), c(0.04, 0.042, 0.044, 0.045))
  result <- yc_zspread(price = 95, coupon_rate = 0.04, maturity = 5,
                        curve = curve, frequency = 2)
  expect_equal(result$model_price, 95, tolerance = 1e-6)
})

test_that("yc_zspread works with annual frequency", {
  curve <- yc_curve(c(1, 2, 5, 10), c(0.04, 0.042, 0.044, 0.045))
  result <- yc_zspread(price = 98, coupon_rate = 0.04, maturity = 5,
                        curve = curve, frequency = 1)
  expect_true(result$zspread > 0)
  expect_equal(result$model_price, 98, tolerance = 1e-6)
})

test_that("yc_zspread accepts list-style curve", {
  curve <- list(maturities = c(1, 2, 5, 10), rates = c(0.04, 0.042, 0.044, 0.045))
  result <- yc_zspread(price = 95, coupon_rate = 0.04, maturity = 5,
                        curve = curve, frequency = 2)
  expect_true(result$zspread > 0)
})

test_that("yc_zspread rejects bad inputs", {
  curve <- yc_curve(c(1, 2, 5), c(0.04, 0.042, 0.044))
  expect_error(yc_zspread(price = -1, coupon_rate = 0.04, maturity = 5,
                           curve = curve), "positive")
  expect_error(yc_zspread(price = 95, coupon_rate = 0.04, maturity = 5,
                           curve = curve, frequency = 4), "frequency")
})

# Key rate duration tests

test_that("yc_key_rate_duration returns correct structure", {
  curve <- yc_curve(c(1, 2, 5, 10, 30), c(0.03, 0.035, 0.04, 0.042, 0.045))
  result <- yc_key_rate_duration(coupon_rate = 0.04, maturity = 10,
                                  curve = curve)
  expect_true(is.data.frame(result))
  expect_true(all(c("tenor", "key_rate_duration") %in% names(result)))
  expect_equal(nrow(result), 5)
})

test_that("yc_key_rate_duration: sum approximates total modified duration", {
  curve <- yc_curve(c(1, 2, 5, 10, 30), c(0.04, 0.042, 0.044, 0.045, 0.046))

  # Price the bond at the curve
  maturity <- 10
  freq <- 2
  step <- 1 / freq
  cf_times <- seq(step, maturity, by = step)
  coupon <- 100 * 0.05 / freq
  cfs <- rep(coupon, length(cf_times))
  cfs[length(cfs)] <- cfs[length(cfs)] + 100
  z <- approx(c(1, 2, 5, 10, 30), c(0.04, 0.042, 0.044, 0.045, 0.046),
              xout = cf_times, rule = 2)$y
  price <- sum(cfs * (1 + z)^(-cf_times))

  krd <- yc_key_rate_duration(coupon_rate = 0.05, maturity = 10,
                               curve = curve, key_rates = c(1, 2, 5, 10, 30))

  # Compute total modified duration via parallel shift
  shift <- 0.0001
  z_up <- z + shift
  price_up <- sum(cfs * (1 + z_up)^(-cf_times))
  total_mod_dur <- -(price_up - price) / (shift * price)

  sum_krd <- sum(krd$key_rate_duration)

  # Sum of KRDs should be close to total modified duration
  expect_equal(sum_krd, total_mod_dur, tolerance = 0.5)
})

test_that("yc_key_rate_duration: all values non-negative for standard bond", {
  curve <- yc_curve(c(1, 2, 5, 10, 30), c(0.04, 0.042, 0.044, 0.045, 0.046))
  result <- yc_key_rate_duration(coupon_rate = 0.05, maturity = 10,
                                  curve = curve)
  # For a standard bond, all KRDs should be non-negative
  expect_true(all(result$key_rate_duration >= -0.01))
})

test_that("yc_key_rate_duration: short bond has KRD concentrated at short end", {
  curve <- yc_curve(c(1, 2, 5, 10, 30), c(0.04, 0.042, 0.044, 0.045, 0.046))
  result <- yc_key_rate_duration(coupon_rate = 0.05, maturity = 2,
                                  curve = curve)
  # For a 2Y bond, KRDs at 10Y and 30Y should be near zero
  expect_true(abs(result$key_rate_duration[result$tenor == 10]) < 0.01)
  expect_true(abs(result$key_rate_duration[result$tenor == 30]) < 0.01)
})

test_that("yc_key_rate_duration rejects bad inputs", {
  curve <- yc_curve(c(1, 2, 5), c(0.04, 0.042, 0.044))
  expect_error(yc_key_rate_duration(coupon_rate = 0.04,
                                     maturity = -5, curve = curve), "positive")
  expect_error(yc_key_rate_duration(coupon_rate = 0.04,
                                     maturity = 5, curve = curve,
                                     frequency = 4), "frequency")
})
