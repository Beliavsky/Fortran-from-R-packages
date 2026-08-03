# ML Tabular Quickstart Example
# Demonstrates ml_backtest() with simple linear models
# Run with: run_example("ml_tabular_quickstart.R")
# Author: Alberto Pallotta
# Note: No external ML libraries required - uses stats::lm.fit()
# ML Tabular Quickstart (easy, offline)
suppressPackageStartupMessages(library(PortfolioTesteR))
suppressPackageStartupMessages(library(data.table))

# 1) Data
data(sample_prices_weekly)
data(sample_prices_daily)

# 2) Features (lagged to avoid look-ahead)
mom   <- panel_lag(calc_momentum(sample_prices_weekly, 12L), 1L)
vol_d <- calc_rolling_volatility(sample_prices_daily, 20L)
vol_w <- panel_lag(align_to_timeframe(vol_d, sample_prices_weekly$Date, "forward_fill"), 1L)

# 3) Labels (future 4-week log return)
Y <- make_labels(sample_prices_weekly, horizon = 4L, type = "log")

# Quick leakage sanity
stopifnot(validate_no_leakage(mom,  Y, verbose = FALSE))
stopifnot(validate_no_leakage(vol_w, Y, verbose = FALSE))

# 4) Simple linear model wrappers
fit_lm  <- function(X, y) { Xc <- cbind(1, as.matrix(X)); list(coef = stats::lm.fit(Xc, y)$coefficients) }
pred_lm <- function(m, X)  as.numeric(cbind(1, as.matrix(X)) %*% m$coef)

# 5) One-call backtest wrapper
set.seed(42)
res <- ml_backtest(
  features_list = list(mom = mom, vol = vol_w),
  labels        = Y,
  fit_fn        = fit_lm,
  predict_fn    = pred_lm,
  schedule      = list(is = 104L, oos = 4L, step = 4L),
  transform     = "zscore",
  selection     = list(top_k = 10L),
  weighting     = list(method = "softmax", temperature = 12),
  caps          = list(max_per_symbol = 0.08),
  prices        = sample_prices_weekly,
  initial_capital = 1e5,
  name          = "Tabular: mom + vol → softmax"
)

print(res$backtest)
if (interactive()) plot(res$backtest, type = "performance")
