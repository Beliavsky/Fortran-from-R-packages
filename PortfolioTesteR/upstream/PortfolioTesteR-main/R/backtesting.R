# backtesting.R
# Core backtesting engine for portfolio strategies
# Designed for extensibility and simplicity

#library(data.table)

###############################################################################
# MAIN BACKTESTING FUNCTION
###############################################################################
#' Run Portfolio Backtest
#'
#' @description
#' Main backtesting engine that simulates portfolio performance over time.
#' Handles position tracking, transaction costs, and performance calculation.
#'
#' @param prices Price data (data.frame with Date column)
#' @param weights Weight matrix from weighting functions
#' @param initial_capital Starting capital (default: 100000)
#' @param name Strategy name for reporting
#' @param verbose Print progress messages (default: FALSE)
#' @param stop_loss Optional stop loss percentage as decimal
#' @param stop_monitoring_prices Optional daily prices for stop monitoring
#'
#' @return backtest_result object with performance metrics
#' @export
#' @examples
#' data("sample_prices_weekly")
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
#' selected <- filter_top_n(momentum, n = 10)
#' weights <- weight_equally(selected)
#' result <- run_backtest(sample_prices_weekly, weights)
run_backtest <- function(prices, weights, initial_capital = 100000,
                         name = "Strategy", verbose = FALSE,
                         stop_loss = NULL, stop_monitoring_prices = NULL) {
  # Run a backtest given prices and pre-calculated weights
  #
  # This is the core backtesting engine that simulates portfolio performance over time.
  # It handles position management, cash tracking, transaction recording, and now includes
  # optional stop loss functionality for risk management.
  #
  # Args:
  #   prices (DataFrame): Price data with Date column and symbol columns
  #                       Should match the frequency of your strategy (daily, weekly, etc.)
  #   weights (DataFrame): Weight allocations with same structure as prices
  #                        Values should sum to 1 for fully invested, 0 for cash
  #   initial_capital (numeric): Starting portfolio value (default: 100000)
  #   name (string): Strategy name for identification (default: "Strategy")
  #   verbose (logical): Print progress information (default: FALSE)
  #   stop_loss (numeric): Optional stop loss as decimal (e.g., 0.15 for 15% stop loss)
  #                        Set to NULL to disable (default: NULL)
  #   stop_monitoring_prices (DataFrame): Daily prices for stop loss monitoring
  #                                      Required if stop_loss is set
  #                                      Must include all symbols in prices
  #                                      Must cover the full backtest period
  #
  # Returns:
  #   list: Backtest results object with class "backtest_result" containing:
  #     - positions: DataFrame of shares held each period
  #     - portfolio_values: Vector of portfolio values over time
  #     - transactions: DataFrame of all trades executed
  #     - returns: Period returns
  #     - stop_loss_trades: Count of positions exited via stop loss
  #     - Comprehensive performance metrics and metadata
  #
  # Stop Loss Behavior:
  #   - Monitors positions daily between rebalance dates
  #   - Triggers when position falls below entry price by stop_loss percentage
  #   - Exits to cash immediately upon trigger
  #   - Prevents re-entry until next rebalance period
  #   - Records stop loss transactions separately for analysis
  #
  # Note:
  #   - Backward compatible: works exactly as before when stop_loss = NULL
  #   - Stop losses are checked on all dates in stop_monitoring_prices between rebalances
  #   - Entry price is set when position is established or increased
  # Validate inputs
  validate_backtest_inputs(prices, weights, initial_capital)

  # NEW: Validate stop loss inputs
  if (!is.null(stop_loss)) {
    if (is.null(stop_monitoring_prices)) {
      stop("stop_monitoring_prices required when stop_loss is set")
    }
    if (stop_loss <= 0 || stop_loss >= 1) {
      stop("stop_loss must be between 0 and 1")
    }
    if (verbose) {
      cat(sprintf("  Stop loss enabled: %.1f%%\n", stop_loss * 100))
    }
  }

  # Convert to data.table for efficiency
  setDT(prices)
  setDT(weights)
  if (!is.null(stop_monitoring_prices)) {
    setDT(stop_monitoring_prices)
  }

  # Get symbol columns
  symbol_cols <- setdiff(names(prices), "Date")
  n_periods <- nrow(prices)

  # AUTO-DETECT WARMUP PERIOD (unchanged)
  weight_sums <- rowSums(weights[, symbol_cols, with = FALSE], na.rm = TRUE)
  first_trade_idx <- which(weight_sums > 0)[1]

  if (is.na(first_trade_idx)) {
    warning("No trades found in strategy!")
    first_trade_idx <- 1
  }

  warmup_periods <- first_trade_idx - 1

  if (verbose) {
    cat(sprintf("Running backtest '%s' with %d symbols over %d periods\n",
                name, length(symbol_cols), n_periods))
    if (warmup_periods > 0) {
      cat(sprintf("  Detected %d warmup periods (no trading)\n", warmup_periods))
    }
  }

  # Initialize result storage
  positions <- data.table(Date = prices$Date)
  portfolio_values <- numeric(n_periods)
  cash <- numeric(n_periods)
  transactions <- list()

  # NEW: Track entry prices and stop status
  entry_prices <- rep(NA_real_, length(symbol_cols))
  names(entry_prices) <- symbol_cols
  position_stopped <- rep(FALSE, length(symbol_cols))
  names(position_stopped) <- symbol_cols

  # Initialize portfolio
  current_shares <- rep(0, length(symbol_cols))
  names(current_shares) <- symbol_cols
  current_cash <- initial_capital

  # Main backtest loop
  for (i in 1:n_periods) {
    date <- prices$Date[i]
    current_prices <- unlist(prices[i, symbol_cols, with = FALSE])
    target_weights <- unlist(weights[i, symbol_cols, with = FALSE])

    # NEW: Apply stop loss if enabled
    if (!is.null(stop_loss) && i > 1) {
      # Find date range for monitoring
      if (i == 1) {
        prev_date <- date
      } else {
        prev_date <- prices$Date[i-1]
      }

      # Get daily prices between rebalances
      monitor_dates <- stop_monitoring_prices$Date[
        stop_monitoring_prices$Date > prev_date &
          stop_monitoring_prices$Date <= date
      ]

      # Check each monitoring day
      for (monitor_date in monitor_dates) {
        monitor_idx <- which(stop_monitoring_prices$Date == monitor_date)
        monitor_prices <- unlist(stop_monitoring_prices[monitor_idx, symbol_cols, with = FALSE])

        # Check stops for held positions
        for (s in seq_along(symbol_cols)) {
          symbol <- symbol_cols[s]

          # Skip if no position or already stopped
          if (current_shares[s] <= 0 || position_stopped[s]) next

          # Check if stop triggered
          if (!is.na(entry_prices[s]) && !is.na(monitor_prices[s])) {
            price_change <- (monitor_prices[s] - entry_prices[s]) / entry_prices[s]

            if (price_change <= -stop_loss) {
              # STOP TRIGGERED!
              if (verbose) {
                cat(sprintf("  Stop loss triggered for %s on %s (%.1f%% loss)\n",
                            symbol, monitor_date, price_change * 100))
              }

              # Execute stop loss sale
              shares_to_sell <- current_shares[s]
              sale_proceeds <- shares_to_sell * monitor_prices[s]

              # Record stop loss transaction
              transactions[[length(transactions) + 1]] <- data.table(
                Date = monitor_date,
                Symbol = symbol,
                Shares = -shares_to_sell,
                Price = monitor_prices[s],
                Value = -sale_proceeds,
                Type = "StopLoss"
              )

              # Update position
              current_shares[s] <- 0
              current_cash <- current_cash + sale_proceeds
              position_stopped[s] <- TRUE

              # Clear target weight to prevent re-entry this period
              target_weights[s] <- 0
            }
          }
        }
      }
    }

    # Calculate current portfolio value
    position_values <- current_shares * current_prices
    position_values[is.na(position_values)] <- 0
    total_position_value <- sum(position_values)
    current_value <- total_position_value + current_cash

    # Calculate target positions (accounting for stopped positions)
    target_values <- current_value * target_weights
    target_shares <- floor(target_values / current_prices)
    target_shares[is.na(target_shares)] <- 0

    # Calculate trades needed
    trades_needed <- target_shares - current_shares

    # Execute trades if any
    if (any(trades_needed != 0)) {

      # STEP 1: PROCESS ALL SELLS FIRST
      sell_value <- 0

      for (s in seq_along(symbol_cols)) {
        symbol <- symbol_cols[s]
        trade_amount <- trades_needed[s]

        if (trade_amount < 0 && current_shares[s] > 0) {  # This is a sell
          shares_to_sell <- abs(trade_amount)
          sell_proceeds <- shares_to_sell * current_prices[s]
          sell_value <- sell_value + sell_proceeds

          # Record transaction
          transactions[[length(transactions) + 1]] <- data.table(
            Date = date,
            Symbol = symbol,
            Shares = -shares_to_sell,
            Price = current_prices[s],
            Value = -sell_proceeds,
            Type = "Rebalance"
          )

          # Update position
          current_shares[s] <- current_shares[s] - shares_to_sell

          # Clear entry price if fully exited
          if (current_shares[s] == 0) {
            entry_prices[s] <- NA_real_
            position_stopped[s] <- FALSE
          }
        }
      }

      # Update cash after all sells
      current_cash <- current_cash + sell_value

      # STEP 2: PROCESS ALL BUYS
      buy_cost <- 0

      for (s in seq_along(symbol_cols)) {
        symbol <- symbol_cols[s]
        trade_amount <- trades_needed[s]

        if (trade_amount > 0) {  # This is a buy
          shares_to_buy <- trade_amount
          cost <- shares_to_buy * current_prices[s]

          # Check if we have enough cash
          if (cost > current_cash + 0.01) {
            next
          }

          buy_cost <- buy_cost + cost

          # Record transaction
          transactions[[length(transactions) + 1]] <- data.table(
            Date = date,
            Symbol = symbol,
            Shares = shares_to_buy,
            Price = current_prices[s],
            Value = cost,
            Type = "Rebalance"
          )

          # Update position and entry price
          current_shares[s] <- current_shares[s] + shares_to_buy
          entry_prices[s] <- current_prices[s]  # NEW: Track entry price
          position_stopped[s] <- FALSE

          # Update cash immediately
          current_cash <- current_cash - cost
        }
      }
    }

    # Record state
    positions[i, (symbol_cols) := as.list(current_shares)]
    cash[i] <- current_cash
    portfolio_values[i] <- current_value

    if (verbose && i %% 10 == 0) {
      cat(sprintf("  Processed %d/%d periods (%.1f%%)\r",
                  i, n_periods, 100 * i / n_periods))
    }
  }

  if (verbose) cat("\n")

  # Combine transactions
  if (length(transactions) > 0) {
    transactions_df <- rbindlist(transactions)
  } else {
    transactions_df <- data.table(
      Date = as.Date(character()),
      Symbol = character(),
      Shares = numeric(),
      Price = numeric(),
      Value = numeric(),
      Type = character()  # NEW: Track transaction type
    )
  }

  # Calculate returns
  returns <- c(0, diff(portfolio_values) / head(portfolio_values, -1))

  # Calculate metrics from first trade
  if (warmup_periods > 0 && first_trade_idx <= n_periods) {
    active_start_value <- portfolio_values[first_trade_idx]
    active_end_value <- tail(portfolio_values, 1)
    active_return <- (active_end_value / active_start_value) - 1
    active_returns <- returns[first_trade_idx:length(returns)]
  } else {
    active_return <- (tail(portfolio_values, 1) / initial_capital) - 1
    active_returns <- returns
  }

  # Create result object
  result <- list(
    # Metadata
    name = name,
    initial_capital = initial_capital,
    start_date = min(prices$Date),
    end_date = max(prices$Date),
    n_periods = n_periods,
    symbols = symbol_cols,

    # NEW: Stop loss info
    stop_loss = stop_loss,
    stop_loss_trades = if (!is.null(stop_loss)) {
      sum(transactions_df$Type == "StopLoss", na.rm = TRUE)
    } else 0,

    # Warmup information
    warmup_periods = warmup_periods,
    active_periods = n_periods - warmup_periods,
    first_trade_date = if(warmup_periods > 0) prices$Date[first_trade_idx] else prices$Date[1],

    # Core data
    dates = prices$Date,
    positions = positions,
    cash = cash,
    portfolio_values = portfolio_values,
    returns = returns,
    transactions = transactions_df,

    # Input data
    prices = prices,
    weights = weights,

    # Summary statistics
    final_value = tail(portfolio_values, 1),
    total_return = active_return,
    total_return_with_warmup = (tail(portfolio_values, 1) / initial_capital) - 1,
    n_transactions = nrow(transactions_df),

    # Store for other calculations
    first_trade_idx = first_trade_idx,
    active_returns = active_returns
  )

  # Add S3 class for custom methods
  class(result) <- c("backtest_result", "list")

  return(result)
}



###############################################################################
# VALIDATION FUNCTIONS
###############################################################################

validate_backtest_inputs <- function(prices, weights, initial_capital) {
  # Validate inputs for backtesting

  # Check data types
  if (!is.data.frame(prices)) {
    stop("prices must be a data.frame or data.table")
  }

  if (!is.data.frame(weights)) {
    stop("weights must be a data.frame or data.table")
  }

  # Check initial capital
  if (!is.numeric(initial_capital) || length(initial_capital) != 1 || initial_capital <= 0) {
    stop("initial_capital must be a positive number")
  }

  # Check structure match
  if (!identical(dim(prices), dim(weights))) {
    stop("prices and weights must have the same dimensions")
  }

  if (!identical(names(prices), names(weights))) {
    stop("prices and weights must have the same column names")
  }

  if (!all(prices$Date == weights$Date)) {
    stop("prices and weights must have the same dates")
  }

  # Check Date column exists
  if (!"Date" %in% names(prices)) {
    stop("prices must have a Date column")
  }

  # Check for valid weights
  symbol_cols <- setdiff(names(weights), "Date")
  weight_sums <- rowSums(weights[, symbol_cols, with = FALSE], na.rm = TRUE)

  if (any(weight_sums > 1.01)) {  # Allow small numerical errors
    problem_dates <- weights$Date[weight_sums > 1.01]
    warning(sprintf("Weights sum to > 1 on %d dates. First date: %s",
                    length(problem_dates), problem_dates[1]))
  }

  return(TRUE)
}

###############################################################################
# S3 METHODS FOR BACKTEST RESULTS
###############################################################################


#' Print Backtest Results
#'
#' @description
#' S3 print method for backtest results. Shows key performance metrics.
#'
#' @param x backtest_result object
#' @param ... Additional arguments (unused)
#'
#' @return Invisible copy of x
#' @export
#' @examples
#' data("sample_prices_weekly")
#' mom <- calc_momentum(sample_prices_weekly, lookback = 12)
#' sel <- filter_top_n(mom, n = 10)
#' W   <- weight_equally(sel)
#' res <- run_backtest(sample_prices_weekly, W)
#' print(res)
print.backtest_result <- function(x, ...) {
  cat("Backtest Result: ", x$name, "\n")
  cat("=====================================\n")

  # Show warmup information if present
  if (x$warmup_periods > 0) {
    cat(sprintf("Warmup Period: %d observations (no trading)\n", x$warmup_periods))
    cat(sprintf("Active Period: %s to %s (%d observations)\n",
                x$first_trade_date, x$end_date, x$active_periods))
  } else {
    cat(sprintf("Period: %s to %s (%d observations)\n",
                x$start_date, x$end_date, x$n_periods))
  }

  cat(sprintf("Initial Capital: $%s\n", format(x$initial_capital, big.mark = ",")))
  cat(sprintf("Final Value: $%s\n", format(round(x$final_value), big.mark = ",")))

  # Show both returns if there's warmup
  if (x$warmup_periods > 0) {
    cat(sprintf("Total Return (active period): %.2f%%\n", x$total_return * 100))
    cat(sprintf("Total Return (full period): %.2f%%\n", x$total_return_with_warmup * 100))
  } else {
    cat(sprintf("Total Return: %.2f%%\n", x$total_return * 100))
  }

  cat(sprintf("Transactions: %d\n", x$n_transactions))

  # Basic performance metrics - calculated on active period
  if (x$active_periods > 0) {
    # Use active period for calculations
    years <- x$active_periods / 52  # Assuming weekly data
    annualized_return <- (1 + x$total_return)^(1/years) - 1

    # Use active returns for volatility
    volatility <- sd(x$active_returns, na.rm = TRUE) * sqrt(52)
    sharpe <- annualized_return / volatility
  } else {
    # Fallback to original calculation
    annualized_return <- calculate_annualized_return(x)
    volatility <- sd(x$returns, na.rm = TRUE) * sqrt(252/7)
    sharpe <- annualized_return / volatility
  }

  cat(sprintf("\nAnnualized Return: %.2f%%\n", annualized_return * 100))
  cat(sprintf("Annualized Volatility: %.2f%%\n", volatility * 100))
  cat(sprintf("Sharpe Ratio: %.2f\n", sharpe))

  # Drawdown info
  dd <- calculate_drawdowns(x$portfolio_values)
  cat(sprintf("Max Drawdown: %.2f%%\n", min(dd) * 100))

  invisible(x)
}

#' Summary method for backtest results
#'
#' @param object A backtest_result object
#' @param ... Additional arguments (unused)
#'
#' @return Invisible copy of the object
#' @export
summary.backtest_result <- function(object, ...) {
  # Detailed summary of backtest results

  x <- object
  cat("\nDetailed Summary: ", x$name, "\n")
  cat("=====================================\n")

  # Position summary
  positions_held <- rowSums(x$positions[, -"Date"] > 0)
  cat("\nPosition Statistics:\n")
  cat(sprintf("  Average positions held: %.1f\n", mean(positions_held)))
  cat(sprintf("  Max positions held: %d\n", max(positions_held)))
  cat(sprintf("  Periods fully invested: %d (%.1f%%)\n",
              sum(x$cash < x$initial_capital * 0.01),
              100 * sum(x$cash < x$initial_capital * 0.01) / x$n_periods))

  # Transaction summary
  if (nrow(x$transactions) > 0) {
    cat("\nTransaction Summary:\n")
    cat(sprintf("  Total trades: %d\n", nrow(x$transactions)))
    cat(sprintf("  Avg trades per period: %.1f\n",
                nrow(x$transactions) / x$n_periods))

    # Turnover calculation
    buy_value <- sum(abs(x$transactions[Shares > 0]$Value))
    avg_portfolio_value <- mean(x$portfolio_values)
    turnover <- buy_value / avg_portfolio_value / (x$n_periods / 52)  # Annualized
    cat(sprintf("  Annual turnover: %.1f%%\n", turnover * 100))
  }

  # Return distribution
  cat("\nReturn Distribution:\n")
  cat(sprintf("  Mean return: %.3f%%\n", mean(x$returns) * 100))
  cat(sprintf("  Median return: %.3f%%\n", median(x$returns) * 100))
  cat(sprintf("  Best period: %.2f%%\n", max(x$returns) * 100))
  cat(sprintf("  Worst period: %.2f%%\n", min(x$returns) * 100))
  cat(sprintf("  Positive periods: %d (%.1f%%)\n",
              sum(x$returns > 0),
              100 * sum(x$returns > 0) / length(x$returns)))

  invisible(x)
}



#' Plot Backtest Results
#'
#' @description
#' S3 plot method for visualizing backtest performance.
#'
#' @param x backtest_result object
#' @param type Plot type: "performance", "drawdown", "weights", or "all"
#' @param ... Additional plotting parameters
#'
#' @return NULL (creates plot)
#' @export
#' @examples
#' data("sample_prices_weekly")
#' mom <- calc_momentum(sample_prices_weekly, lookback = 12)
#' sel <- filter_top_n(mom, n = 10)
#' W   <- weight_equally(sel)
#' res <- run_backtest(sample_prices_weekly, W)
#' if (interactive()) plot(res, type = "performance")
plot.backtest_result <- function(x, type = "performance", ...) {
  # Save and restore only mfrow
  old_mfrow <- par("mfrow")
  on.exit(par(mfrow = old_mfrow))

  # Plot backtest results
  if (type == "performance") {
    # Equity curve
    par(mfrow = c(2, 1))

    # Portfolio value
    plot(x$dates, x$portfolio_values, type = "l",
         main = paste(x$name, "- Portfolio Value"),
         xlab = "Date", ylab = "Portfolio Value ($)",
         col = "darkblue", lwd = 2)
    grid()

    # Drawdown
    dd <- calculate_drawdowns(x$portfolio_values)
    plot(x$dates, dd * 100, type = "l",
         main = "Drawdown %",
         xlab = "Date", ylab = "Drawdown (%)",
         col = "darkred", lwd = 2)
    grid()
    abline(h = 0, lty = 2)

  } else if (type == "positions") {
    # Position counts over time
    positions_held <- rowSums(x$positions[, -"Date"] > 0)
    plot(x$dates, positions_held, type = "l",
         main = paste(x$name, "- Number of Positions"),
         xlab = "Date", ylab = "Number of Holdings",
         col = "darkgreen", lwd = 2)
    grid()

  } else if (type == "weights") {
    # Weight distribution over time
    symbol_cols <- setdiff(names(x$weights), "Date")
    # Get top 5 by average weight
    avg_weights <- colMeans(x$weights[, symbol_cols, with = FALSE], na.rm = TRUE)
    top_symbols <- names(sort(avg_weights, decreasing = TRUE))[1:min(5, length(avg_weights))]
    # Plot weights
    plot(x$dates, rep(0, length(x$dates)), type = "n",
         ylim = c(0, 1),
         main = paste(x$name, "- Portfolio Weights"),
         xlab = "Date", ylab = "Weight")
    colors <- c("red", "blue", "green", "orange", "purple")
    for (i in seq_along(top_symbols)) {
      lines(x$dates, x$weights[[top_symbols[i]]],
            col = colors[i], lwd = 2)
    }
    legend("topright", legend = top_symbols,
           col = colors[1:length(top_symbols)], lwd = 2)
    grid()
  }
}

###############################################################################
# UTILITY FUNCTIONS
###############################################################################

#' Calculate Portfolio Drawdowns
#'
#' @description
#' Calculates drawdown series from portfolio values. Drawdown is the
#' percentage decline from the previous peak.
#'
#' @param values Numeric vector of portfolio values
#'
#' @return Numeric vector of drawdown percentages (negative values)
#' @keywords internal
calculate_drawdowns <- function(values) {
  # FIXED: Handle zero/NA in cummax
  cummax_values <- cummax(values)

  # Avoid division by zero
  drawdowns <- safe_divide(values - cummax_values, cummax_values)

  # If cummax is 0 (shouldn't happen with positive portfolio values), set to 0
  drawdowns[cummax_values == 0] <- 0

  return(drawdowns)
}


#' Calculate Annualized Return
#'
#' @description
#' Converts total return to annualized return based on time period.
#'
#' @param result Backtest result object
#'
#' @return Annualized return as decimal (0.1 = 10%)
#' @keywords internal
calculate_annualized_return <- function(result) {
  # Calculate annualized return
  years <- as.numeric(result$end_date - result$start_date) / 365.25
  annualized <- (result$final_value / result$initial_capital)^(1/years) - 1
  return(annualized)
}


#' Calculate Comprehensive Backtest Metrics
#'
#' @description
#' Computes performance metrics including Sharpe ratio, maximum drawdown,
#' win rate, and other statistics from backtest results.
#'
#' @param result Backtest result object from run_backtest()
#'
#' @return List containing performance metrics
#' @export
#' @examples
#' # Create a backtest result to use
#' data(sample_prices_weekly)
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
#' selected <- filter_top_n(momentum, n = 10)
#' weights <- weight_equally(selected)
#' result <- run_backtest(sample_prices_weekly, weights)
#'
#' # Calculate metrics
#' metrics <- backtest_metrics(result)
#' print(metrics$sharpe_ratio)
backtest_metrics <- function(result) {
  # Calculate comprehensive metrics

  # Use active returns if available
  returns <- if (!is.null(result$active_returns)) result$active_returns else result$returns
  values <- result$portfolio_values

  # Basic metrics - use total_return which is now based on active period
  total_return <- result$total_return

  # Calculate annualized return based on active periods
  if (result$warmup_periods > 0) {
    years <- result$active_periods / periods_per_year
    annualized_return <- (1 + total_return)^(1/years) - 1
  } else {
    years <- as.numeric(result$end_date - result$start_date) / 365.25
    annualized_return <- (result$final_value / result$initial_capital)^(1/years) - 1
  }

  # Risk metrics - calculated on active returns
  volatility <- sd(returns, na.rm = TRUE)
  downside_returns <- returns[returns < 0]
  downside_vol <- sd(downside_returns, na.rm = TRUE)

  # Annualize based on frequency (assume weekly for now)
  # Detect frequency and annualize correctly
  frequency <- get_data_frequency(result$dates)
  periods_per_year <- switch(frequency,
                             "daily" = 252,
                             "weekly" = 52,
                             "monthly" = 12,
                             "quarterly" = 4,
                             52  # default to weekly if unknown
  )
  ann_vol <- volatility * sqrt(periods_per_year)
  ann_downside_vol <- downside_vol * sqrt(periods_per_year)

  # Risk-adjusted returns
  sharpe <- annualized_return / ann_vol
  sortino <- annualized_return / ann_downside_vol

  # Drawdown metrics - calculated on full period (this is correct)
  dd <- calculate_drawdowns(values)
  max_dd <- min(dd)

  # Win/loss statistics - on active returns
  winning_periods <- sum(returns > 0)
  total_periods <- length(returns)
  win_rate <- winning_periods / total_periods

  avg_win <- mean(returns[returns > 0])
  avg_loss <- mean(returns[returns < 0])

  metrics <- list(
    total_return = total_return,
    annualized_return = annualized_return,
    volatility = ann_vol,
    sharpe_ratio = sharpe,
    sortino_ratio = sortino,
    max_drawdown = max_dd,
    win_rate = win_rate,
    avg_win = avg_win,
    avg_loss = avg_loss,
    total_trades = result$n_transactions,
    warmup_periods = result$warmup_periods,
    active_periods = result$active_periods
  )

  return(metrics)
}
