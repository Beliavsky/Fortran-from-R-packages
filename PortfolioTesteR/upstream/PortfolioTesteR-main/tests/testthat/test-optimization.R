test_that("run_param_grid returns a well-formed result and enforces guardrails", {
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

  opt <- run_param_grid(
    prices  = sample_prices_weekly,
    grid    = grid,
    builder = builder,
    metric  = NULL,          # default metric_sharpe()
    verbose = FALSE,
    light_mode = TRUE,
    precompute_returns = TRUE
  )

  # structure
  expect_s3_class(opt, "param_grid_result")
  expect_true(is.list(opt$all_results))
  expect_equal(length(opt$all_scores), nrow(do.call(expand.grid, grid)))
  expect_true(is.numeric(opt$optimization_score) || is.na(opt$optimization_score))
  expect_true(is.null(opt$best_params) || is.list(opt$best_params))
  expect_true(is.null(opt$best_idx) || (opt$best_idx >= 1 && opt$best_idx <= length(opt$all_results)))

  # printing should produce output (not be "silent")
  expect_output(print(opt), "Parameter Grid Optimization Result", fixed = TRUE)

  # plotting should not error (render to a null device)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(opt, type = "line",    params = "lookback"))
  expect_no_error(plot(opt, type = "heatmap", params = c("lookback", "n_top")))
})

test_that("run_param_grid rejects cadence/timeframe knobs", {
  skip_on_cran()
  data(sample_prices_weekly)

  builder <- function(prices, params, ...) {
    weight_equally(filter_top_n(calc_momentum(prices, params$lookback), 5))
  }

  # Any forbidden name should error (e.g., "rebalance")
  expect_error(
    run_param_grid(
      prices  = sample_prices_weekly,
      grid    = list(lookback = c(8, 12), rebalance = c(4, 8)),
      builder = builder
    ),
    "forbidden parameter names",
    ignore.case = TRUE
  )
})

