test_that("yc_curve creates valid object", {
  d <- make_simple_yields()
  curve <- yc_curve(d$maturities, d$rates)
  expect_s3_class(curve, "yc_curve")
  expect_equal(curve$method, "observed")
  expect_equal(curve$type, "zero")
  expect_equal(curve$n_obs, 4)
  expect_null(curve$fitted)
  expect_null(curve$residuals)
})

test_that("yc_curve sorts by maturity", {
  curve <- yc_curve(c(10, 2, 30, 5), c(0.04, 0.045, 0.043, 0.042))
  expect_equal(curve$maturities, c(2, 5, 10, 30))
  expect_equal(curve$rates, c(0.045, 0.042, 0.04, 0.043))
})

test_that("yc_curve accepts date", {
  curve <- yc_curve(c(2, 10), c(0.04, 0.05), date = as.Date("2024-01-15"))
  expect_equal(curve$date, as.Date("2024-01-15"))
})

test_that("yc_curve accepts type argument", {
  curve <- yc_curve(c(2, 10), c(0.04, 0.05), type = "par")
  expect_equal(curve$type, "par")
})

test_that("yc_curve rejects non-numeric maturities", {
  expect_error(yc_curve("a", 0.04), "numeric")
})

test_that("yc_curve rejects negative maturities", {
  expect_error(yc_curve(c(-1, 2), c(0.04, 0.05)), "positive")
})

test_that("yc_curve rejects NA in rates", {
  expect_error(yc_curve(c(1, 2), c(NA, 0.05)), "NA")
})

test_that("yc_curve rejects mismatched lengths", {
  expect_error(yc_curve(c(1, 2), c(0.04)), "same length")
})

test_that("yc_curve rejects invalid date", {
  expect_error(yc_curve(c(1, 2), c(0.04, 0.05), date = "2024-01-01"), "Date")
})

test_that("yc_curve print method works", {
  d <- make_simple_yields()
  curve <- yc_curve(d$maturities, d$rates)
  expect_no_error(print(curve))
})
