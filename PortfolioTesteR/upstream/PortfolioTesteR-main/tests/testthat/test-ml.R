# tests/testthat/test-ml-v2.R

test_that("combine_scores and roll_fit_predict run", {
  data(sample_prices_weekly); data(sample_prices_daily)

  mom <- panel_lag(calc_momentum(sample_prices_weekly, 12), 1)
  vol <- panel_lag(
    align_to_timeframe(calc_rolling_volatility(sample_prices_daily, 20),
                       sample_prices_weekly$Date, "forward_fill"),
    1
  )
  Y <- make_labels(sample_prices_weekly, horizon = 4, type = "log")

  fit_lm  <- function(X, y) { Xc <- cbind(1, X); list(coef = stats::lm.fit(Xc, y)$coefficients) }
  pred_lm <- function(m, X) as.numeric(cbind(1, X) %*% m$coef)

  # old:
  # S <- roll_fit_predict_v2(list(mom = mom, vol = vol), Y, fit_lm, pred_lm, 52, 4, 4)

  # new:
  S <- roll_fit_predict(list(mom = mom, vol = vol), Y, fit_lm, pred_lm,
                        is_periods = 52, oos_periods = 4, step = 4)

  expect_s3_class(S, "data.table")
  expect_true("Date" %in% names(S))
  expect_gt(length(setdiff(names(S), "Date")), 0L)
})

