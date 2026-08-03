# example_risk_parity.R
# =====================
# Demonstrates defensive portfolio construction using volatility weighting.
# Lower volatility stocks receive higher allocations.

library(PortfolioTesteR)

# Load sample data
data("sample_prices_weekly")
prices <- sample_prices_weekly

# Calculate momentum for stock selection
momentum <- calc_momentum(prices, lookback = 12)

# Select top 15 stocks by momentum
selected <- filter_top_n(momentum, n = 15)

# Calculate volatility (20-week)
volatility <- calc_rolling_volatility(prices, window = 20)

# Weight inversely by volatility (defensive approach)
# low_vol_preference = TRUE gives more weight to stable stocks
weights <- weight_by_volatility(
  selected_df = selected,
  volatility_df = volatility,
  low_vol_preference = TRUE
)

# Run backtest
result <- run_backtest(
  prices = prices,
  weights = weights,
  initial_capital = 100000,
  name = "Defensive Momentum"
)

print(result)
