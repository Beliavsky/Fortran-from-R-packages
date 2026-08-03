# Optimization surfaces and slices (3D + facets)
library(PortfolioTesteR)
set.seed(1)
data(sample_prices_weekly)

builder_opt3 <- function(prices, params, ...) {
  sig <- calc_momentum(prices, params$lookback)
  if (!is.null(params$thresh)) {
    # optional pre-filter on signal level before top-n
    sig <- filter_above(sig, params$thresh)
  }
  sel <- filter_top_n(sig, params$n_top)
  weight_equally(sel)
}

# small grid but enough for 3D/slices
grid_opt3 <- list(
  lookback = c(8, 12, 16, 20),
  n_top    = c(5, 10, 15),
  thresh   = c(0.00, 0.05, 0.10)  # adjust if your momentum scale differs
)

opt3 <- run_param_grid(
  prices  = sample_prices_weekly,
  grid    = grid_opt3,
  builder = builder_opt3,
  metric  = NULL,
  name_prefix = "MOM3",
  verbose = FALSE,
  light_mode = TRUE,
  precompute_returns = TRUE
)

# 3D surface over two params
plot(opt3, type = "surface", params = c("lookback", "n_top"))

# 2D heatmaps faceted by third param
plot(opt3, type = "slices",
     params = c("lookback", "n_top", "thresh"),
     facet_values = c(0.00, 0.05, 0.10), ncol = 3)

# 3D surfaces faceted by third param
plot(opt3, type = "surface_slices",
     params = c("lookback", "n_top", "thresh"),
     facet_values = c(0.00, 0.05, 0.10), theta = 30, phi = 25)

# Let the method auto-pick a layout (optional)
plot(opt3, type = "auto")
