test_that("yc_fit dispatches to nelson_siegel", {
  d <- make_us_yields()
  fit <- yc_fit(d$maturities, d$rates, method = "nelson_siegel")
  expect_equal(fit$method, "nelson_siegel")
})

test_that("yc_fit dispatches to svensson", {
  d <- make_us_yields()
  fit <- yc_fit(d$maturities, d$rates, method = "svensson")
  expect_equal(fit$method, "svensson")
})

test_that("yc_fit dispatches to cubic_spline", {
  d <- make_us_yields()
  fit <- yc_fit(d$maturities, d$rates, method = "cubic_spline")
  expect_equal(fit$method, "cubic_spline")
})

test_that("yc_fit default is nelson_siegel", {
  d <- make_simple_yields()
  fit <- yc_fit(d$maturities, d$rates)
  expect_equal(fit$method, "nelson_siegel")
})
