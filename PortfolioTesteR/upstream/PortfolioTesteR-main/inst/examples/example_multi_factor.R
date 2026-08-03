# example_multi_factor.R
# ======================
# Demonstrates combining multiple signals and filters.
# Uses momentum, volatility, and trend filters together.

library(PortfolioTesteR)

# Load sample data
data("sample_prices_weekly")
prices <- sample_prices_weekly

# Calculate indicators
momentum <- calc_momentum(prices, lookback = 12)
volatility <- calc_rolling_volatility(prices, window = 20)
ma50 <- calc_moving_average(prices, window = 50)

# Filter 1: Positive momentum
positive_momentum <- filter_above(momentum, value = 0)

# Filter 2: Below median volatility
median_vol <- calc_cross_sectional_percentile(volatility, percentile = 50)
low_volatility <- filter_below(volatility, median_vol)

# Filter 3: Price above MA50
above_trend <- filter_above(calc_distance(prices, ma50), value = 0)

# Combine all filters (all conditions must be true)
all_conditions <- combine_filters(
  list(positive_momentum, low_volatility, above_trend),
  op = "and"
)

# Select top 10 by momentum from qualified stocks
selected <- filter_top_n_where(
  signal_df = momentum,
  n = 10,
  condition_df = all_conditions
)

# Equal weight the final selection
weights <- weight_equally(selected)

# Run backtest
result <- run_backtest(
  prices = prices,
  weights = weights,
  initial_capital = 100000,
  name = "Multi-Factor Strategy"
)

print(result)
