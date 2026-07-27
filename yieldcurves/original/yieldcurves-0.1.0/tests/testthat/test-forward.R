test_that("yc_forward computes instantaneous forwards", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  fwd <- yc_forward(fit)
  expect_equal(nrow(fwd), length(d$maturities))
  expect_true(all(c("maturity", "forward_rate") %in% names(fwd)))
  # Forward rates should be in a reasonable range
  expect_true(all(fwd$forward_rate > 0 & fwd$forward_rate < 0.10))
})

test_that("yc_forward computes forward-forward rates", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  fwd <- yc_forward(fit, maturities = c(1, 2, 5), horizon = 1)
  expect_equal(nrow(fwd), 3)
})

test_that("yc_forward at custom maturities", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)
  fwd <- yc_forward(fit, maturities = c(1, 5, 10))
  expect_equal(fwd$maturity, c(1, 5, 10))
})

test_that("yc_forward NS analytical matches numerical", {
  d <- make_us_yields()
  fit <- yc_nelson_siegel(d$maturities, d$rates)

  # Analytical forward at specific maturities
  fwd_analytical <- yc_forward(fit, maturities = c(2, 5, 10))

  # Numerical forward via finite difference on spot * maturity
  m <- c(2, 5, 10)
  dm <- 1e-6
  r <- yc_predict(fit, m)$rate
  r_plus <- yc_predict(fit, m + dm)$rate
  fwd_numerical <- (r_plus * (m + dm) - r * m) / dm

  expect_equal(fwd_analytical$forward_rate, fwd_numerical, tolerance = 1e-4)
})

test_that("yc_forward works with observed curve", {
  d <- make_simple_yields()
  curve <- yc_curve(d$maturities, d$rates)
  fwd <- yc_forward(curve)
  expect_equal(nrow(fwd), length(d$maturities))
})

test_that("yc_forward Svensson analytical matches numerical", {
  d <- make_us_yields()
  fit <- yc_svensson(d$maturities, d$rates)

  fwd_analytical <- yc_forward(fit, maturities = c(2, 5, 10))

  m <- c(2, 5, 10)
  dm <- 1e-6
  r <- yc_predict(fit, m)$rate
  r_plus <- yc_predict(fit, m + dm)$rate
  fwd_numerical <- (r_plus * (m + dm) - r * m) / dm

  expect_equal(fwd_analytical$forward_rate, fwd_numerical, tolerance = 1e-4)
})

test_that("yc_forward rejects non-curve", {
  expect_error(yc_forward(list()), "yc_curve")
})
