# Basic parameter grid optimization (line plot)
library(PortfolioTesteR)
set.seed(1)
data(sample_prices_weekly)

builder_opt <- function(prices, params, ...) {
  sel <- filter_top_n(calc_momentum(prices, params$lookback), params$n_top)
  weight_equally(sel)
}

grid_opt <- list(lookback = c(8, 12, 26), n_top = c(5, 10))

opt <- run_param_grid(
  prices  = sample_prices_weekly,
  grid    = grid_opt,
  builder = builder_opt,
  metric  = NULL,                 # defaults to metric_sharpe()
  name_prefix = "MOM",
  verbose = FALSE,
  light_mode = TRUE,
  precompute_returns = TRUE
)

print(opt)
plot(opt, type = "line", params = "lookback")
