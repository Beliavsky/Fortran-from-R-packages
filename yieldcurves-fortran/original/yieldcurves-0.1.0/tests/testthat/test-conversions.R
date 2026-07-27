test_that("yc_par_to_zero returns correct structure", {
  result <- yc_par_to_zero(c(1, 2, 3), c(0.04, 0.042, 0.043))
  expect_true(is.data.frame(result))
  expect_true(all(c("maturity", "zero_rate") %in% names(result)))
  expect_equal(nrow(result), 3)
})

test_that("yc_par_to_zero: 1-year par = 1-year zero", {
  result <- yc_par_to_zero(c(1, 2, 3), c(0.04, 0.042, 0.043))
  expect_equal(result$zero_rate[1], 0.04)
})

test_that("yc_par_to_zero: flat par curve gives flat zero curve", {
  result <- yc_par_to_zero(c(1, 2, 3, 5), rep(0.05, 4))
  # For a flat par curve, zero rates should be very close to par
  expect_equal(result$zero_rate, rep(0.05, 4), tolerance = 0.001)
})

test_that("yc_zero_to_par returns correct structure", {
  result <- yc_zero_to_par(c(1, 2, 3), c(0.04, 0.042, 0.043))
  expect_true(is.data.frame(result))
  expect_true(all(c("maturity", "par_rate") %in% names(result)))
})

test_that("yc_zero_to_par: 1-year zero = 1-year par", {
  result <- yc_zero_to_par(c(1, 2, 3), c(0.04, 0.042, 0.043))
  expect_equal(result$par_rate[1], 0.04)
})

test_that("yc_zero_to_par: flat zero curve gives flat par curve", {
  result <- yc_zero_to_par(c(1, 2, 3, 5), rep(0.05, 4))
  expect_equal(result$par_rate, rep(0.05, 4), tolerance = 0.001)
})

test_that("par_to_zero and zero_to_par roundtrip", {
  maturities <- c(1, 2, 3, 5, 10)
  original_par <- c(0.040, 0.042, 0.043, 0.044, 0.045)

  zeros <- yc_par_to_zero(maturities, original_par)
  back_to_par <- yc_zero_to_par(maturities, zeros$zero_rate)

  expect_equal(back_to_par$par_rate, original_par, tolerance = 0.001)
})

test_that("yc_par_to_zero hand-computed 2Y zero rate", {
  # Par rates: 1Y = 4%, 2Y = 4.2%
  # 1Y zero = 4%
  # 2Y: 1 = 0.042/(1.04) + 1.042/(1+z2)^2
  # PV coupon = 0.042/1.04 = 0.04038462
  # (1+z2)^2 = 1.042 / (1 - 0.04038462) = 1.042/0.9596154 = 1.08597
  # z2 = sqrt(1.08597) - 1 = 0.04209
  result <- yc_par_to_zero(c(1, 2), c(0.04, 0.042))
  # Verify: PV = 0.042/1.04 + 1.042/(1+z2)^2 = 1
  z2 <- result$zero_rate[2]
  pv <- 0.042 / 1.04 + 1.042 / (1 + z2)^2
  expect_equal(pv, 1, tolerance = 1e-8)
})

test_that("yc_par_to_zero warns for non-integer maturities", {
  expect_warning(
    yc_par_to_zero(c(1, 2.5), c(0.04, 0.042)),
    "Non-integer"
  )
})

test_that("yc_par_to_zero rejects bad inputs", {
  expect_error(yc_par_to_zero("a", 0.04), "numeric")
  expect_error(yc_par_to_zero(c(1, 2), c(0.04)), "same length")
})
