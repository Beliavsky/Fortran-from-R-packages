# Demo: StochRSI Acceleration + Inverse-Vol Risk Parity
# ------------------------------------------------------
# What this demonstrates
# - Momentum-of-momentum via StochRSI "acceleration" (StochRSI vs its MA)
# - Selection uses a high-StochRSI gate + top-N by acceleration
# - Allocation by inverse-volatility risk parity using DAILY prices
# - Backtest on the weekly grid (bundled datasets only; runs fast)

library(PortfolioTesteR)
set.seed(1)

# --- Data (bundled) -------------------------------------------------------
data(sample_prices_weekly)
data(sample_prices_daily)

# Exclude broad ETFs from stock-selection universe
symbols_all   <- setdiff(names(sample_prices_weekly), "Date")
stock_symbols <- setdiff(symbols_all, c("SPY", "TLT"))

weekly_stocks <- sample_prices_weekly[, c("Date", stock_symbols), with = FALSE]
daily_stocks  <- sample_prices_daily[,  c("Date", stock_symbols), with = FALSE]

# --- StochRSI "acceleration" signal (weekly) ------------------------------
stochrsi    <- calc_stochrsi(weekly_stocks, length = 14)   # in [0, 1]
stochrsi_ma <- calc_moving_average(stochrsi, window = 5)   # short MA of StochRSI
accel       <- calc_distance(stochrsi, stochrsi_ma)        # positive = rising

# Gate to "high StochRSI" zone, then take top-12 by acceleration
high_zone <- filter_above(stochrsi, value = 0.80)
sel      <- filter_top_n_where(
  signal_df     = accel,
  n             = 12,
  condition_df  = high_zone,
  min_qualified = 8,
  ascending     = FALSE
)

# --- Allocation: inverse-volatility risk parity (DAILY prices) ------------
w_ivol <- weight_by_risk_parity(
  selected_df      = sel,
  prices_df        = daily_stocks,
  method           = "inverse_vol",
  lookback_periods = 126,  # ~6 months; small for speed
  min_periods      = 60
)

# --- Backtest on the weekly grid ------------------------------------------
res <- run_backtest(
  prices          = weekly_stocks,
  weights         = w_ivol,
  initial_capital = 100000,
  name            = "StochRSI Accel + InvVol RP"
)

print(res)
summary(res)
# plot(res, type = "performance")  # uncomment when running interactively
