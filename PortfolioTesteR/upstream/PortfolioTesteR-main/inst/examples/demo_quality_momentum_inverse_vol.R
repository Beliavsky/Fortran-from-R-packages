# Demo: Quality-Momentum + Inverse-Vol Risk Parity with Mild Regime
# -----------------------------------------------------------------
# - Quality: RSI in [40, 60] band
# - Selection: from qualified names, pick top-12 momentum; cap to 10 positions
# - Allocation: inverse-volatility risk parity (DAILY prices)
# - Regime: SPY above 20W MA = risk-on; else keep 30% exposure (partial)
# - Backtest on weekly grid; no external I/O

library(PortfolioTesteR)
set.seed(1)

# --- Data (bundled) -------------------------------------------------------
data(sample_prices_weekly)
data(sample_prices_daily)

symbols_all  <- setdiff(names(sample_prices_weekly), "Date")
weekly_spy   <- sample_prices_weekly[, c("Date", "SPY"), with = FALSE]
weekly_trading <- sample_prices_weekly    # use full universe for selection
daily_trading  <- sample_prices_daily

# --- Quality (RSI band) + Momentum (weekly) -------------------------------
rsi    <- calc_rsi(weekly_trading, period = 14)
q_band <- filter_between(rsi, lower = 40, upper = 60)

mom12  <- calc_momentum(weekly_trading, lookback = 12)

# Select top-12 momentum FROM the qualified set; require at least 8 qualified
sel_qm <- filter_top_n_where(
  signal_df     = mom12,
  n             = 12,
  condition_df  = q_band,
  min_qualified = 8,
  ascending     = FALSE
)

# Cap final breadth to 10 positions, ranked by momentum
sel_qm <- limit_positions(
  selection_df   = sel_qm,
  max_positions  = 10,
  ranking_signal = mom12
)

# --- Mild regime (partial) mask via SPY 20W MA ----------------------------
ma20_spy <- calc_moving_average(weekly_spy, window = 20)
dist_ma  <- calc_distance(weekly_spy, ma20_spy)
risk_on  <- dist_ma$SPY > 0

sel_reg <- apply_regime(
  selection_df     = sel_qm,
  regime_condition = risk_on,
  partial_weight   = 0.3      # keep 30% in risk-off instead of going flat
)

# --- Inverse-volatility risk parity (DAILY prices) ------------------------
w_iv <- weight_by_risk_parity(
  selected_df      = sel_reg,
  prices_df        = daily_trading,
  method           = "inverse_vol",
  lookback_periods = 126,
  min_periods      = 60
)

# --- Backtest on weekly grid ----------------------------------------------
res <- run_backtest(
  prices          = weekly_trading,
  weights         = w_iv,
  initial_capital = 100000,
  name            = "Quality-MOM + Inverse-Vol (mild regime)"
)

print(res)
summary(res)
# plot(res, type = "performance")  # uncomment when running interactively



