# Demo: Regime-Filtered Momentum with HRP Allocation (weekly IS, daily risk)
# -------------------------------------------------------------------------
# What it shows
# - SPY trend regime (above/below 20-week MA) gating
# - Cross-sectional momentum selection on the weekly grid
# - HRP (Hierarchical Risk Parity) allocation using DAILY prices
# - Backtest on the weekly grid (no external data; runs in seconds)

library(PortfolioTesteR)
set.seed(1)

# --- Data (bundled) -------------------------------------------------------
data(sample_prices_weekly)
data(sample_prices_daily)

# Use SPY for regime; exclude it from the trading universe
symbols_all <- setdiff(names(sample_prices_weekly), "Date")
trading_symbols <- setdiff(symbols_all, "SPY")

weekly_trading <- sample_prices_weekly[, c("Date", trading_symbols), with = FALSE]
daily_trading  <- sample_prices_daily[,  c("Date", trading_symbols), with = FALSE]
weekly_spy     <- sample_prices_weekly[, c("Date", "SPY"), with = FALSE]

# --- Signals & selection (weekly) ----------------------------------------
mom12 <- calc_momentum(weekly_trading, lookback = 12)
sel10 <- filter_top_n(mom12, n = 10)

# --- Simple trend regime from SPY (weekly) --------------------------------
ma20_spy <- calc_moving_average(weekly_spy, window = 20)
dist_ma  <- calc_distance(weekly_spy, ma20_spy)      # >0 means price above MA
regime   <- dist_ma$SPY > 0                          # TRUE = “risk on”

# Apply regime: fully exit in bad regime (partial_weight = 0)
sel_reg <- apply_regime(
  selection_df    = sel10,
  regime_condition = regime,
  partial_weight   = 0
)

# --- HRP allocation (risk model from DAILY prices) ------------------------
w_hrp <- weight_by_hrp(
  selected_df     = sel_reg,
  prices_df       = daily_trading,
  lookback_periods = 126,
  min_periods      = 60
)

# --- Backtest on the weekly grid -----------------------------------------
res <- run_backtest(
  prices = weekly_trading,
  weights = w_hrp,
  initial_capital = 100000,
  name = "Regime-Filtered MOM + HRP"
)

print(res)
summary(res)
# Uncomment to see a chart when running interactively:
#plot(res, type = "performance")
