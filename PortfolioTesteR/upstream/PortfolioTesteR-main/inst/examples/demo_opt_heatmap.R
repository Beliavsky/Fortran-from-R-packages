# Parameter grid optimization (2D heatmap)
library(PortfolioTesteR)
set.seed(1)
data(sample_prices_weekly)

builder_opt <- function(prices, params, ...) {
  sel <- filter_top_n(calc_momentum(prices, params$lookback), params$n_top)
  weight_equally(sel)
}

# Small, fast grid
grid_opt <- list(lookback = c(8, 12, 16, 20), n_top = c(5, 10, 15))

opt <- run_param_grid(
  prices  = sample_prices_weekly,
  grid    = grid_opt,
  builder = builder_opt,
  metric  = NULL,
  name_prefix = "MOM",
  verbose = FALSE,
  light_mode = TRUE,
  precompute_returns = TRUE
)

print(opt)
plot(opt, type = "heatmap", params = c("lookback", "n_top"))
