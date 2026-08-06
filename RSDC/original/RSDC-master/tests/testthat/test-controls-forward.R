# Regression test: control$NP / control$steptol must reach the global search
# for every estimator, including "const" (they were hardcoded there before
# 1.7-0's audit fix).

test_that("const forwards control$NP and control$steptol to the theta search", {
  skip_on_cran()
  real_search <- RSDC:::.rsdc_theta_search
  seen <- NULL
  testthat::local_mocked_bindings(
    .rsdc_theta_search = function(...) {
      seen <<- list(...)$con
      real_search(...)
    },
    .package = "RSDC"
  )
  set.seed(1)
  y <- matrix(rnorm(160), ncol = 2)
  fit <- rsdc_estimate(method = "const", residuals = y,
                       control = list(itermax = 15, NP = 20, steptol = 5,
                                      compute_se = FALSE))
  expect_s3_class(fit, "rsdc_fit")
  expect_identical(seen$NP, 20)
  expect_identical(seen$steptol, 5)
  expect_identical(seen$itermax, 15)
})
