# example_momentum_basic.R
# ========================
# Demonstrates the core PortfolioTesteR workflow:
# Calculate indicator → Filter stocks → Weight portfolio → Run backtest
#
# This example selects the 10 strongest momentum stocks and weights them equally.
# Perfect starting point for learning the library's pattern.

library(PortfolioTesteR)

# Load sample data included with the package
data("sample_prices_weekly")
prices <- sample_prices_weekly

# Step 1: Calculate momentum (12-week returns)
# The lookback period determines the timeframe for measuring performance
momentum <- calc_momentum(prices, lookback = 12)

# Step 2: Select top 10 stocks by momentum
# This creates a binary matrix: 1 = selected, 0 = not selected
selected <- filter_top_n(momentum, n = 10)

# Step 3: Weight the selected stocks
# Equal weighting gives each selected stock 10% of the portfolio
weights <- weight_equally(selected)

# Step 4: Run the backtest
# This simulates buying and rebalancing according to the weights
result <- run_backtest(
  prices = prices,
  weights = weights,
  initial_capital = 100000,
  name = "Simple Momentum"
)

# View results
print(result)

# Optional: Calculate performance metrics
metrics <- backtest_metrics(result)
print(sprintf("Annual Return: %.1f%%", metrics$annualized_return * 100))
print(sprintf("Sharpe Ratio: %.2f", metrics$sharpe_ratio))
print(sprintf("Max Drawdown: %.1f%%", metrics$max_drawdown * 100))
