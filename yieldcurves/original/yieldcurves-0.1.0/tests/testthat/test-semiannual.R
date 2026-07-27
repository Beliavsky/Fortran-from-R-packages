test_that("par_to_zero with frequency=1 works identically to original", {
  maturities <- c(1, 2, 3, 5, 10)
  par_rates <- c(0.040, 0.042, 0.043, 0.044, 0.045)
  result <- yc_par_to_zero(maturities, par_rates, frequency = 1)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 5)
  # 1Y zero = par

  expect_equal(result$zero_rate[1], 0.04)
})

test_that("zero_to_par with frequency=1 works identically to original", {
  maturities <- c(1, 2, 3, 5, 10)
  zero_rates <- c(0.040, 0.042, 0.043, 0.044, 0.045)
  result <- yc_zero_to_par(maturities, zero_rates, frequency = 1)
  expect_true(is.data.frame(result))
  expect_equal(result$par_rate[1], 0.04)
})

test_that("semi-annual bootstrap: known 2Y 5% bond", {
  # 2-year bond paying 5% semi-annually
  # Cash flows at 0.5, 1.0, 1.5, 2.0
  # If par rates are flat at 5%, zero rates should also be ~5%
  result <- yc_par_to_zero(c(0.5, 1, 1.5, 2), rep(0.05, 4), frequency = 2)
  expect_equal(result$zero_rate, rep(0.05, 4), tolerance = 0.001)
})

test_that("semi-annual par_to_zero and zero_to_par roundtrip", {
  maturities <- c(0.5, 1, 1.5, 2, 3, 5)
  par_rates <- c(0.040, 0.042, 0.043, 0.044, 0.046, 0.048)

  zeros <- yc_par_to_zero(maturities, par_rates, frequency = 2)
  back <- yc_zero_to_par(maturities, zeros$zero_rate, frequency = 2)

  expect_equal(back$par_rate, par_rates, tolerance = 0.001)
})

test_that("semi-annual bootstrap: 0.5Y par = 0.5Y zero", {
  result <- yc_par_to_zero(c(0.5, 1, 2), c(0.04, 0.042, 0.043), frequency = 2)
  expect_equal(result$zero_rate[1], 0.04)
})

test_that("frequency parameter rejects invalid values", {
  expect_error(
    yc_par_to_zero(c(1, 2), c(0.04, 0.042), frequency = 4),
    "frequency"
  )
  expect_error(
    yc_zero_to_par(c(1, 2), c(0.04, 0.042), frequency = 3),
    "frequency"
  )
})

test_that("annual roundtrip still works with explicit frequency=1", {
  maturities <- c(1, 2, 3, 5, 10)
  original_par <- c(0.040, 0.042, 0.043, 0.044, 0.045)

  zeros <- yc_par_to_zero(maturities, original_par, frequency = 1)
  back_to_par <- yc_zero_to_par(maturities, zeros$zero_rate, frequency = 1)

  expect_equal(back_to_par$par_rate, original_par, tolerance = 0.001)
})
