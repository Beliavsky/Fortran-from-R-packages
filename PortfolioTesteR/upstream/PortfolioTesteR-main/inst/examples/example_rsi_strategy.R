# example_rsi_strategy.R
# ======================
# Demonstrates mean reversion using RSI with signal-based weighting.
# Selects stocks in moderate RSI range and weights by signal strength.

library(PortfolioTesteR)

# Load sample data
data("sample_prices_weekly")
prices <- sample_prices_weekly

# Calculate RSI (14-period standard)
rsi <- calc_rsi(prices, period = 14)

# Select stocks with RSI between 40-60 (avoiding extremes)
selected <- filter_between(rsi, lower = 40, upper = 60)

# Weight by RSI signal strength
# Stocks with higher RSI (within selected range) get more weight
weights <- weight_by_signal(selected, rsi)

# Run backtest
result <- run_backtest(
  prices = prices,
  weights = weights,
  initial_capital = 100000,
  name = "RSI Mean Reversion"
)

print(result)
