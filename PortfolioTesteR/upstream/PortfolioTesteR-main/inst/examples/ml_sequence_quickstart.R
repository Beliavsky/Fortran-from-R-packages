# ML Sequence Quickstart
# Demonstrates ml_backtest_seq() with sequence features
# Run with: run_example("ml_sequence_quickstart.R")
# No external ML libraries required - uses stats::lm.fit()
# ML Sequence Quickstart (medium, offline)
suppressPackageStartupMessages(library(PortfolioTesteR))
suppressPackageStartupMessages(library(data.table))

# 1) Data
data(sample_prices_weekly)
data(sample_prices_daily)

# 2) Base tabular features (lag to avoid leakage)
mom   <- panel_lag(calc_momentum(sample_prices_weekly, 12L), 1L)
vol_d <- calc_rolling_volatility(sample_prices_daily, 20L)
vol_w <- panel_lag(align_to_timeframe(vol_d, sample_prices_weekly$Date, "forward_fill"), 1L)

# 3) Labels for a 4-week horizon
Y <- make_labels(sample_prices_weekly, horizon = 4L, type = "log")
stopifnot(validate_no_leakage(mom,  Y, verbose = FALSE))
stopifnot(validate_no_leakage(vol_w, Y, verbose = FALSE))

# 4) Sequence model (simple linear fit on flattened sequences)
fit_seq_lm  <- function(X, y) { Xc <- cbind(1, as.matrix(X)); list(coef = stats::lm.fit(Xc, y)$coefficients) }
pred_seq_lm <- function(m, X)  as.numeric(cbind(1, as.matrix(X)) %*% m$coef)

# 5) Backtest with sequences
set.seed(123)
res_seq <- ml_backtest_seq(
  features_list   = list(mom = mom, vol = vol_w),
  labels          = Y,
  steps           = 26L, horizon = 4L,
  fit_fn          = fit_seq_lm,
  predict_fn      = pred_seq_lm,
  schedule        = list(is = 104L, oos = 4L, step = 4L),
  group           = "pooled",
  normalize       = "zscore",
  selection       = list(top_k = 15L),
  weighting       = list(method = "softmax", temperature = 12),
  caps            = list(max_per_symbol = 0.08),
  prices          = sample_prices_weekly,
  initial_capital = 1e5,
  name            = "SEQ: 26x2 → softmax"
)

print(res_seq$backtest)
if (interactive()) plot(res_seq$backtest, type = "performance")

# # Optional: tiny purged/embargoed CV inside IS (kept small)
# is_start <- 1L; is_end <- 104L
# cv <- cv_tune_seq(
#   features_list = list(mom = mom, vol = vol_w),
#   labels        = Y,
#   is_start      = is_start,
#   is_end        = is_end,
#   steps_grid    = c(12L, 26L),
#   horizon       = 4L,
#   fit_fn        = fit_seq_lm,
#   predict_fn    = pred_seq_lm,
#   k             = 2L,
#   purge         = NULL,
#   embargo       = NULL,
#   group         = "pooled",
#   max_train_samples = 2000L,
#   max_val_samples   = 1000L,
#   na_action     = "omit",
#   metric        = "spearman"
# )
# print(cv)

