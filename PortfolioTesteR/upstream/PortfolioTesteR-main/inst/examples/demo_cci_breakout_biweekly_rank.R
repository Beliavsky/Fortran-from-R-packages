# Demo: CCI Breakout with Bi-Weekly Cadence + Rank Weights
# --------------------------------------------------------
# What this demonstrates
# - Breakout filter using CCI > +100 on weekly prices
# - From breakers, select top-10 by CCI; weight by rank (exponential)
# - Bi-weekly RE-TUNING of weights; carry forward to weekly dates
# - Backtest on weekly grid (bundled datasets; fast)

library(PortfolioTesteR)
set.seed(1)

# --- Data (bundled) -------------------------------------------------------
data(sample_prices_weekly)

# Exclude broad ETFs from stock-selection
symbols_all   <- setdiff(names(sample_prices_weekly), "Date")
stock_symbols <- setdiff(symbols_all, c("SPY", "TLT"))

weekly_stocks <- sample_prices_weekly[, c("Date", stock_symbols), with = FALSE]

# --- CCI breakout + rank weighting (weekly) -------------------------------
cci <- calc_cci(weekly_stocks, period = 20)

breakout <- filter_above(cci, value = 100)   # classic CCI breakout threshold

sel <- filter_top_n_where(
  signal_df     = cci,
  n             = 10,
  condition_df  = breakout,
  min_qualified = 6,
  ascending     = FALSE
)

w_rank <- weight_by_rank(
  selected_df = sel,
  signal_df   = cci,
  method      = "exponential",
  ascending   = FALSE
)

# --- Bi-weekly cadence: recompute every 2 weeks, hold between ----------------
# Keep every second rebalancing row, then align forward to weekly dates
# (weights persist between decision dates)
w_bi <- w_rank[seq(1, nrow(w_rank), by = 2), ]
w_aligned <- align_to_timeframe(
  high_freq_data = w_bi,                  # sparser (bi-weekly) decisions
  low_freq_dates = weekly_stocks$Date,    # weekly backtest grid
  method         = "forward_fill"
)

# --- Backtest on the weekly grid --------------------------------------------
res <- run_backtest(
  prices          = weekly_stocks,
  weights         = w_aligned,
  initial_capital = 100000,
  name            = "CCI Breakout (bi-weekly, exp-rank)"
)

print(res)
summary(res)
# plot(res, type = "performance")  # uncomment when running interactively
