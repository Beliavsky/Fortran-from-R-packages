# Walk-forward full diagnostics (windows, drawdown, distributions)
library(PortfolioTesteR)
set.seed(1)
data(sample_prices_weekly)

builder_wf <- function(prices, params, ...) {
  weight_equally(
    filter_top_n(
      calc_momentum(prices, params$lookback),
      params$n_top
    )
  )
}

grid_wf <- list(lookback = c(8, 12, 26), n_top = c(5, 10))

wf <- run_walk_forward(
  prices      = sample_prices_weekly,
  grid        = grid_wf,
  builder     = builder_wf,
  metric      = NULL,      # defaults to metric_sharpe()
  is_periods  = 52,
  oos_periods = 26,
  step        = 26,
  verbose     = FALSE,
  light_mode  = TRUE,
  precompute_all = TRUE,
  builder_args = list()
)


plot(wf, type = "windows",     metric = "OOS_Score")
plot(wf, type = "drawdown")
plot(wf, type = "distributions")
