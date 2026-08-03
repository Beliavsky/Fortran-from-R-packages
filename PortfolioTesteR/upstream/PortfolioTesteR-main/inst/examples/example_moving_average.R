# example_moving_average.R
# ========================
# Demonstrates trend following using moving averages.
# Buy stocks trading above their moving average, avoid those below.

library(PortfolioTesteR)

# Load sample data
data("sample_prices_weekly")
prices <- sample_prices_weekly

# Calculate 20-week moving average
ma20 <- calc_moving_average(prices, window = 20)

# Calculate distance from MA (positive = above, negative = below)
distance <- calc_distance(prices, ma20)

# Select stocks trading above their MA
above_ma <- filter_above(distance, value = 0)

# From those above MA, select top 10 by momentum
momentum <- calc_momentum(prices, lookback = 8)
selected <- filter_top_n_where(
  signal_df = momentum,
  n = 10,
  condition_df = above_ma
)

# Equal weight the selected stocks
weights <- weight_equally(selected)

# Run backtest
result <- run_backtest(
  prices = prices,
  weights = weights,
  initial_capital = 100000,
  name = "MA Trend Following"
)

print(result)
