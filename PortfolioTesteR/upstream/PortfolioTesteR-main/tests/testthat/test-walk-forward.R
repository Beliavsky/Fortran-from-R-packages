test_that("run_walk_forward creates consistent windows and stitched OOS equity", {
  skip_on_cran()

  data(sample_prices_weekly)

  builder <- function(prices, params, ...) {
    weight_equally(
      filter_top_n(
        calc_momentum(prices, params$lookback),
        params$n_top
      )
    )
  }

  grid <- list(lookback = c(8, 12), n_top = c(5, 10))

  wf <- run_walk_forward(
    prices      = sample_prices_weekly,
    grid        = grid,
    builder     = builder,
    metric      = NULL,      # default metric_sharpe()
    is_periods  = 52,
    oos_periods = 26,
    step        = 26,
    verbose     = FALSE,
    light_mode  = TRUE,
    precompute_all = TRUE,
    builder_args = list()
  )

  # structure
  expect_s3_class(wf, "wf_optimization_result")
  expect_true(is.list(wf$windows))
  expect_true(length(wf$windows) >= 1L)

  # NOTE: the object uses 'chosen_params' (not 'best_params')
  expect_true(is.list(wf$chosen_params))
  expect_equal(length(wf$chosen_params), length(wf$windows))

  # stitched OOS equity is a data.frame with Date/Value
  expect_true(is.data.frame(wf$oos_stitched))
  expect_true(all(c("Date", "Value") %in% names(wf$oos_stitched)))
  expect_true(nrow(wf$oos_stitched) >= 1L)

  # optimization summary has expected fields
  expect_true(is.data.frame(wf$optimization_summary))
  expect_true(all(c("Window","IS_Score","OOS_Score","OOS_Return") %in% names(wf$optimization_summary)))

  # plot smoke tests (to null device)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(wf, type = "parameters"))
  expect_no_error(plot(wf, type = "is_oos", metric = "OOS_Score"))
  expect_no_error(plot(wf, type = "equity"))
})
