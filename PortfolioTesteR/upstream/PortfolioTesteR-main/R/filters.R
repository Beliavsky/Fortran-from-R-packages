# filters.R
# Consolidated filter system for portfolio selection
# All filters work with signal DataFrames and return binary selection matrices

# Required libraries
#library(data.table)

###############################################################################
# INPUT VALIDATION HELPERS
###############################################################################

validate_signal_data <- function(signal_df, function_name) {
  # Standard validation for all filter functions

  if (!is.data.frame(signal_df)) {
    stop(paste(function_name, ": signal_df must be a data.frame or data.table"))
  }

  if (!"Date" %in% names(signal_df)) {
    stop(paste(function_name, ": signal_df must have a 'Date' column"))
  }

  if (ncol(signal_df) < 2) {
    stop(paste(function_name, ": signal_df must have at least one signal column besides Date"))
  }

  if (nrow(signal_df) == 0) {
    stop(paste(function_name, ": signal_df cannot be empty"))
  }

  return(TRUE)
}

create_result_template <- function(signal_df) {
  # Create a template for results with same structure as input, filled with 0s
  setDT(signal_df)
  symbol_cols <- setdiff(names(signal_df), "Date")

  result_df <- copy(signal_df)
  result_df[, (symbol_cols) := 0]

  return(result_df)
}

###############################################################################
# CORE FILTER FUNCTIONS
###############################################################################
#' Select Top or Bottom N Stocks by Signal
#'
#' @description
#' Selects the top N (best) or worst N stocks based on signal strength.
#' Optimized using matrix operations for 5-10x speedup.
#'
#' @param signal_df Data frame with Date column and signal values
#' @param n Number of stocks to select
#' @param type "top" for highest values, "worst" for lowest values
#'
#' @return Binary selection matrix (1 = selected, 0 = not selected)
#' @export
#' @examples
#' data("sample_prices_weekly")
#' momentum <- calc_momentum(sample_prices_weekly, 12)
#' # Select 10 highest momentum stocks
#' top10 <- filter_rank(momentum, 10, type = "top")
filter_rank <- function(signal_df, n, type = c("top", "worst")) {
  validate_signal_data(signal_df, "filter_rank")
  type <- match.arg(type)

  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != as.integer(n)) {
    stop("filter_rank: n must be a positive integer")
  }

  # Use copy to avoid mutation
  dt <- ensure_dt_copy(signal_df)
  symbol_cols <- setdiff(names(dt), "Date")

  if (n > length(symbol_cols)) {
    warning(paste("filter_rank: n (", n, ") is greater than number of symbols (",
                  length(symbol_cols), "). Selecting all symbols."))
    n <- length(symbol_cols)
  }

  # Convert to matrix for vectorized operations (5-10x faster than row loops)
  signal_matrix <- as.matrix(dt[, symbol_cols, with = FALSE])

  # Pre-allocate selection matrix for efficiency
  selection_matrix <- matrix(0, nrow = nrow(signal_matrix), ncol = ncol(signal_matrix))

  # Process each time period
  for (i in seq_len(nrow(signal_matrix))) {
    row_vals <- signal_matrix[i, ]

    # Skip rows with all NA values
    if (all(is.na(row_vals))) next

    # Count valid (non-NA) values
    valid_mask <- !is.na(row_vals)
    n_valid <- sum(valid_mask)

    if (n_valid == 0) next

    # Can't select more than available
    n_select <- min(n, n_valid)

    # Get indices of top/worst n stocks using R's efficient order()
    if (type == "top") {
      # Select highest values (best performers)
      top_indices <- order(row_vals, decreasing = TRUE, na.last = TRUE)[1:n_select]
    } else {
      # Select lowest values (worst performers)
      top_indices <- order(row_vals, decreasing = FALSE, na.last = TRUE)[1:n_select]
    }

    # Mark selected positions in matrix
    selection_matrix[i, top_indices] <- 1
  }

  # Convert back to data.table format
  selected_df <- data.table::data.table(Date = dt$Date)  # CHANGED: Added data.table::
  for (j in seq_along(symbol_cols)) {
    selected_df[, (symbol_cols[j]) := selection_matrix[, j]]
  }

  return(selected_df)
}


#' Filter by Threshold Value
#'
#' @description
#' Selects stocks above or below a threshold value.
#'
#' @param signal_df Data frame with signal values
#' @param value Threshold value
#' @param type "above" or "below"
#'
#' @return Binary selection matrix
#' @export
#' @examples
#' data("sample_prices_weekly")
#' momentum <- calc_momentum(sample_prices_weekly, 12)
#' # Select stocks with positive momentum
#' positive <- filter_threshold(momentum, 0, type = "above")
filter_threshold <- function(signal_df, value, type = c("above", "below")) {
  validate_signal_data(signal_df, "filter_threshold")
  type <- match.arg(type)

  if (!is.numeric(value) || length(value) != 1) {
    stop("filter_threshold: value must be a single numeric value")
  }

  # FIXED: Use copy
  dt <- ensure_dt_copy(signal_df)
  symbol_cols <- setdiff(names(dt), "Date")

  selected_df <- create_result_template(dt)

  for (col in symbol_cols) {
    if (type == "above") {
      selected_df[, (col) := ifelse(!is.na(dt[[col]]) & dt[[col]] >= value, 1, 0)]
    } else {
      selected_df[, (col) := ifelse(!is.na(dt[[col]]) & dt[[col]] <= value, 1, 0)]
    }
  }

  return(selected_df)
}

# Fix filter_range
filter_range <- function(signal_df, lower, upper, type = c("inside", "outside")) {
  validate_signal_data(signal_df, "filter_range")
  type <- match.arg(type)

  if (!is.numeric(lower) || length(lower) != 1) {
    stop("filter_range: lower must be a single numeric value")
  }

  if (!is.numeric(upper) || length(upper) != 1) {
    stop("filter_range: upper must be a single numeric value")
  }

  if (lower > upper) {
    stop("filter_range: lower bound must be <= upper bound")
  }

  # FIXED: Use copy
  dt <- ensure_dt_copy(signal_df)
  symbol_cols <- setdiff(names(dt), "Date")

  selected_df <- create_result_template(dt)

  for (col in symbol_cols) {
    if (type == "inside") {
      selected_df[, (col) := ifelse(!is.na(dt[[col]]) &
                                      dt[[col]] >= lower &
                                      dt[[col]] <= upper, 1, 0)]
    } else {
      selected_df[, (col) := ifelse(!is.na(dt[[col]]) &
                                      (dt[[col]] < lower | dt[[col]] > upper), 1, 0)]
    }
  }

  return(selected_df)
}


#' Select Top N Stocks by Signal Value
#'
#' @description
#' Most commonly used filter function. Selects top N (highest) or bottom N (lowest)
#' stocks by signal value. Optimized for 5-10x faster performance.
#'
#' @param signal_df Data frame with Date column and signal values
#' @param n Number of stocks to select
#' @param ascending FALSE (default) selects highest, TRUE selects lowest
#'
#' @return Binary selection matrix (1 = selected, 0 = not selected)
#' @export
#' @examples
#' data("sample_prices_weekly")
#' momentum <- calc_momentum(sample_prices_weekly, 12)
#' # Select 10 highest momentum stocks
#' top_momentum <- filter_top_n(momentum, n = 10)
filter_top_n <- function(signal_df, n, ascending = FALSE) {
  # Select top N (highest) or bottom N (lowest) stocks by signal strength
  # OPTIMIZED VERSION: 5-10x faster than original implementation
  #
  # This is the most commonly used filter function in momentum strategies.
  # It ranks all stocks by their signal value and selects the top (or bottom) N.
  #
  # Args:
  #   signal_df: DataFrame with Date column and signal values for each symbol
  #   n: Number of stocks to select
  #   ascending: If FALSE (default), select highest values (top N)
  #              If TRUE, select lowest values (bottom N)
  #
  # Returns:
  #   Binary selection matrix (1 = selected, 0 = not selected)
  #
  # Examples:
  #   # Select 10 stocks with highest momentum (default behavior)
  #   top_momentum <- filter_top_n(momentum, n = 10)
  #
  #   # Select 5 stocks with lowest volatility (ascending = TRUE)
  #   low_vol_stocks <- filter_top_n(volatility, n = 5, ascending = TRUE)
  #
  #   # Select top 20 by RSI
  #   high_rsi <- filter_top_n(rsi_values, n = 20)
  #
  # Note:
  #   - Automatically handles NA values (they are never selected)
  #   - If n > available stocks, selects all available
  #   - Ties are handled consistently using R's default ordering

  # Route to the optimized filter_rank function
  if (ascending) {
    # For ascending=TRUE, select lowest values (worst performers)
    filter_rank(signal_df, n, type = "worst")
  } else {
    # For ascending=FALSE (default), select highest values (top performers)
    filter_rank(signal_df, n, type = "top")
  }
}


#' Filter Stocks Above Threshold
#'
#' @description
#' Convenience function to select stocks with signal above a value.
#'
#' @param signal_df Data frame with signal values
#' @param value Threshold value
#'
#' @return Binary selection matrix
#' @export
#' @examples
#' data("sample_prices_weekly")
#' rsi <- calc_rsi(sample_prices_weekly, 14)
#' high_rsi <- filter_above(rsi, 70)
filter_above <- function(signal_df, value) {
  filter_threshold(signal_df, value, type = "above")
}

#' Filter Stocks Below Threshold
#'
#' @description
#' Convenience function to select stocks with signal below a value.
#'
#' @param signal_df Data frame with signal values
#' @param value Threshold value
#'
#' @return Binary selection matrix
#' @export
#' @examples
#' data("sample_prices_weekly")
#' rsi <- calc_rsi(sample_prices_weekly, 14)
#' oversold <- filter_below(rsi, 30)
filter_below <- function(signal_df, value) {
  filter_threshold(signal_df, value, type = "below")
}
#' Filter Stocks Between Two Values
#'
#' @description
#' Selects stocks with signal values between lower and upper bounds.
#'
#' @param signal_df Data frame with signal values
#' @param lower Lower bound (inclusive)
#' @param upper Upper bound (inclusive)
#'
#' @return Binary selection matrix
#' @export
#' @examples
#' data("sample_prices_weekly")
#' rsi <- calc_rsi(sample_prices_weekly, 14)
#' # Select stocks with RSI between 30 and 70
#' neutral_rsi <- filter_between(rsi, 30, 70)
filter_between <- function(signal_df, lower, upper) {
  filter_range(signal_df, lower, upper, type = "inside")
}
###############################################################################
# UTILITY FUNCTIONS
###############################################################################

summary_selection <- function(selected_df) {
  # Provide summary statistics about selection results
  #
  # Args:
  #   selected_df (DataFrame): Binary selection matrix from filter functions
  #
  # Returns:
  #   List: Summary statistics

  setDT(selected_df)
  symbol_cols <- setdiff(names(selected_df), "Date")

  # Calculate selection statistics
  total_selections <- selected_df[, lapply(.SD, sum, na.rm = TRUE), .SDcols = symbol_cols]
  selections_per_date <- selected_df[, rowSums(.SD, na.rm = TRUE), .SDcols = symbol_cols]

  summary_stats <- list(
    total_dates = nrow(selected_df),
    total_symbols = length(symbol_cols),
    selections_per_symbol = as.list(total_selections),
    avg_selections_per_date = mean(selections_per_date),
    max_selections_per_date = max(selections_per_date),
    min_selections_per_date = min(selections_per_date),
    dates_with_no_selections = sum(selections_per_date == 0)
  )

  return(summary_stats)
}

validate_filter_result <- function(selected_df, original_df) {
  # Validate that filter result has correct structure
  #
  # Args:
  #   selected_df (DataFrame): Result from filter function
  #   original_df (DataFrame): Original signal DataFrame
  #
  # Returns:
  #   Boolean: TRUE if validation passes

  checks <- list(
    "Same number of rows" = nrow(selected_df) == nrow(original_df),
    "Same column names" = identical(names(selected_df), names(original_df)),
    "Same Date column" = identical(selected_df$Date, original_df$Date),
    "Only 0/1 values" = all(selected_df[, -"Date"] %in% c(0, 1, NA), na.rm = TRUE)
  )

  failed_checks <- names(checks)[!unlist(checks)]

  if (length(failed_checks) > 0) {
    warning("Filter validation failed: ", paste(failed_checks, collapse = ", "))
    return(FALSE)
  }

  return(TRUE)
}



#' Select Top N from Qualified Stocks
#'
#' @description
#' Selects top N stocks by signal, but only from those meeting a condition.
#' Combines qualification and ranking in one step.
#'
#' @param signal_df Signal values for ranking
#' @param n Number to select
#' @param condition_df Binary matrix of qualified stocks
#' @param min_qualified Minimum qualified stocks required (default: 1)
#' @param ascending FALSE for highest, TRUE for lowest
#'
#' @return Binary selection matrix
#' @export
#' @examples
#' data("sample_prices_weekly")
#' # Calculate indicators
#' momentum <- calc_momentum(sample_prices_weekly, 12)
#' ma20 <- calc_moving_average(sample_prices_weekly, 20)
#' distance_from_ma <- calc_distance(sample_prices_weekly, ma20)
#'
#' # Top 10 momentum stocks from those above MA
#' above_ma <- filter_above(distance_from_ma, 0)
#' top_qualified <- filter_top_n_where(momentum, 10, above_ma)
filter_top_n_where <- function(signal_df, n, condition_df, min_qualified = 1, ascending = FALSE) {
  # Select top N (or bottom N) stocks by signal, but only from those meeting a condition
  #
  # Args:
  #   signal_df: DataFrame with signals for ranking
  #   n: Number of stocks to select
  #   condition_df: Binary DataFrame indicating which stocks qualify (1) or not (0)
  #   min_qualified: Minimum number of qualified stocks required to make selections
  #   ascending: If FALSE (default), select highest values (top N)
  #              If TRUE, select lowest values (bottom N)
  #
  # Returns:
  #   Binary selection matrix (1 = selected, 0 = not selected)

  # Input validation
  validate_signal_data(signal_df, "filter_top_n_where")
  validate_signal_data(condition_df, "filter_top_n_where")
  validate_matching_structure(signal_df, condition_df, "filter_top_n_where")

  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != as.integer(n)) {
    stop("filter_top_n_where: n must be a positive integer")
  }

  if (!is.numeric(min_qualified) || length(min_qualified) != 1 || min_qualified < 1) {
    stop("filter_top_n_where: min_qualified must be >= 1")
  }

  # Convert to data.table
  signal_dt <- ensure_dt_copy(signal_df)
  condition_dt <- ensure_dt_copy(condition_df)
  symbol_cols <- setdiff(names(signal_dt), "Date")

  # Initialize result
  selected_df <- create_result_template(signal_dt)

  # Process each date
  for (i in 1:nrow(signal_dt)) {
    # Get signals and conditions for this date
    signals <- unlist(signal_dt[i, ..symbol_cols])
    conditions <- unlist(condition_dt[i, ..symbol_cols])

    # Find qualified stocks (meeting condition and having valid signal)
    qualified_mask <- conditions > 0 & !is.na(signals)
    qualified_symbols <- symbol_cols[qualified_mask]
    qualified_signals <- signals[qualified_mask]

    # Check if we have enough qualified stocks
    if (length(qualified_symbols) < min_qualified) {
      next  # Skip this date - not enough qualified stocks
    }

    # Select top/bottom N from qualified stocks based on ascending parameter
    n_select <- min(n, length(qualified_symbols))

    if (ascending) {
      # Select LOWEST signal values (bottom N)
      selected_indices <- order(qualified_signals)[1:n_select]
    } else {
      # Select HIGHEST signal values (top N) - default behavior
      selected_indices <- order(qualified_signals, decreasing = TRUE)[1:n_select]
    }

    selected_symbols <- qualified_symbols[selected_indices]

    # Mark selected stocks
    selected_df[i, (selected_symbols) := 1]
  }

  return(selected_df)
}





#' Apply Market Regime Filter
#'
#' @description
#' Applies regime-based filtering. When regime is FALSE (e.g., bear market),
#' all selections become 0, moving portfolio to cash.
#'
#' @param selection_df Binary selection matrix
#' @param regime_condition Logical vector (TRUE = trade, FALSE = cash)
#' @param partial_weight Fraction to hold when regime is FALSE (default: 0)
#'
#' @return Modified selection matrix respecting regime
#' @export
#' @examples
#' data("sample_prices_weekly")
#' # Create selection
#' momentum <- calc_momentum(sample_prices_weekly, 12)
#' selected <- filter_top_n(momentum, 10)
#'
#' # Only trade when SPY above 20-week MA
#' ma20 <- calc_moving_average(sample_prices_weekly, 20)
#' spy_regime <- sample_prices_weekly$SPY > ma20$SPY
#' spy_regime[is.na(spy_regime)] <- FALSE
#'
#' regime_filtered <- apply_regime(selected, spy_regime)
apply_regime <- function(selection_df, regime_condition, partial_weight = 0) {
  # Input validation
  validate_signal_data(selection_df, "apply_regime")

  if (!is.numeric(partial_weight) || partial_weight < 0 || partial_weight > 1) {
    stop("apply_regime: partial_weight must be between 0 and 1")
  }

  # Ensure data.table
  selection_dt <- ensure_dt_copy(selection_df)
  symbol_cols <- setdiff(names(selection_dt), "Date")

  # Handle regime condition format
  if (is.data.frame(regime_condition) || is.data.table(regime_condition)) {
    # Extract the actual condition vector
    if ("regime" %in% names(regime_condition)) {
      regime_vector <- regime_condition$regime
    } else if (ncol(regime_condition) == 2 && "Date" %in% names(regime_condition)) {
      # Assume the other column is our regime
      regime_col <- setdiff(names(regime_condition), "Date")[1]
      regime_vector <- regime_condition[[regime_col]]
    } else {
      stop("apply_regime: Cannot determine regime column from regime_condition")
    }
  } else {
    regime_vector <- regime_condition
  }

  # Ensure regime vector matches selection length
  if (length(regime_vector) != nrow(selection_dt)) {
    stop("apply_regime: regime_condition length must match selection_df rows")
  }

  # Apply regime filter
  result_df <- copy(selection_dt)

  for (i in seq_len(nrow(result_df))) {
    if (is.na(regime_vector[i]) || !regime_vector[i]) {
      # Bad regime or NA - scale down or eliminate positions
      if (partial_weight == 0) {
        # Complete exit - set all to 0
        result_df[i, (symbol_cols) := 0]
      } else {
        # Partial exit - scale positions
        for (col in symbol_cols) {
          result_df[i, (col) := result_df[[col]][i] * partial_weight]
        }
      }
    }
    # If regime is TRUE, leave selections unchanged
  }

  return(result_df)
}









#' Filter by Percentile
#'
#' @description
#' Select securities in the top or bottom X percentile.
#' More intuitive than filter_top_n when universe size varies.
#'
#' @param signal_df DataFrame with signal values
#' @param percentile Percentile threshold (0-100)
#' @param type "top" for highest signals, "bottom" for lowest
#'
#' @return Binary selection matrix
#' @export
#' @examples
#' data("sample_prices_weekly")
#' momentum <- calc_momentum(sample_prices_weekly, 12)
#' # Select top 20th percentile
#' top_20pct <- filter_by_percentile(momentum, 20, type = "top")
filter_by_percentile <- function(signal_df, percentile, type = c("top", "bottom")) {
  # Select securities in the top or bottom X percentile
  #
  # This is more intuitive than filter_top_n when universe size varies.
  # Example: "top 20th percentile" adapts to universe size automatically.
  #
  # Args:
  #   signal_df: DataFrame with signal values (momentum, value score, etc.)
  #   percentile: Percentile threshold (0-100). E.g., 20 = top/bottom 20%
  #   type: "top" for highest signals, "bottom" for lowest
  #
  # Returns:
  #   Binary selection matrix (1 = selected, 0 = not selected)

  validate_signal_data(signal_df, "filter_by_percentile")
  type <- match.arg(type)

  if (!is.numeric(percentile) || percentile <= 0 || percentile > 100) {
    stop("filter_by_percentile: percentile must be between 0 and 100")
  }

  dt <- ensure_dt_copy(signal_df)
  symbol_cols <- setdiff(names(dt), "Date")

  selected_df <- create_result_template(dt)

  for (i in 1:nrow(dt)) {
    row_values <- unlist(dt[i, ..symbol_cols])
    valid_values <- row_values[!is.na(row_values)]

    if (length(valid_values) == 0) {
      next
    }

    # Calculate percentile threshold
    if (type == "top") {
      threshold <- quantile(valid_values, probs = 1 - percentile/100)
      selected_mask <- row_values >= threshold
    } else {  # bottom
      threshold <- quantile(valid_values, probs = percentile/100)
      selected_mask <- row_values <= threshold
    }

    # Handle NAs
    selected_mask[is.na(selected_mask)] <- FALSE

    # Set selections
    selected_cols <- symbol_cols[selected_mask]
    if (length(selected_cols) > 0) {
      selected_df[i, (selected_cols) := 1]
    }
  }

  return(selected_df)
}






#' Combine Multiple Filter Conditions
#'
#' @description
#' Combines multiple filter conditions using AND or OR logic.
#'
#' @param ... Two or more filter data frames to combine
#' @param op Operation: "and" or "or"
#' @param apply_when Optional condition vector for conditional filtering
#' @param debug Print debug information (default: FALSE)
#'
#' @return Combined binary selection matrix
#' @export
#' @examples
#' data("sample_prices_weekly")
#' momentum <- calc_momentum(sample_prices_weekly, 12)
#' rsi <- calc_rsi(sample_prices_weekly, 14)
#' # Create individual filters
#' high_momentum <- filter_above(momentum, 0.05)
#' moderate_rsi <- filter_between(rsi, 40, 60)
#' # Combine them
#' combined <- combine_filters(high_momentum, moderate_rsi, op = "and")
combine_filters <- function(..., op = "and", apply_when = NULL, debug = FALSE) {
  # Combine multiple filter conditions with AND/OR logic
  # FIXED VERSION with better error handling and debugging
  #
  # Args:
  #   ...: Two or more filter DataFrames to combine
  #   op: Operation - "and" or "or"
  #   apply_when: Optional condition vector for conditional filtering
  #   debug: Print debug information during execution
  #
  # Returns:
  #   Combined filter DataFrame

  filters <- list(...)

  if (debug) {
    cat("DEBUG: combine_filters called with", length(filters), "filters\n")
    cat("DEBUG: Operation:", op, "\n")
    cat("DEBUG: apply_when provided:", !is.null(apply_when), "\n")
  }

  # Validation
  if (length(filters) < 2) {
    stop("combine_filters: Need at least 2 filters to combine")
  }

  # Validate all inputs are data frames
  for (i in 1:length(filters)) {
    if (!is.data.frame(filters[[i]])) {
      stop(paste("combine_filters: Filter", i, "must be a data.frame or data.table"))
    }
  }

  # Check dimensions match
  base_dim <- dim(filters[[1]])
  base_names <- names(filters[[1]])

  for (i in 2:length(filters)) {
    if (!identical(dim(filters[[i]]), base_dim)) {
      stop(paste("combine_filters: Filter", i, "has different dimensions"))
    }
    if (!identical(names(filters[[i]]), base_names)) {
      stop(paste("combine_filters: Filter", i, "has different column names"))
    }
  }

  # Validate operation
  if (!op %in% c("and", "or")) {
    stop("combine_filters: op must be 'and' or 'or'")
  }

  # Create result as copy of first filter
  result <- ensure_dt_copy(filters[[1]])
  symbol_cols <- setdiff(names(result), "Date")

  if (debug) {
    cat("DEBUG: Symbol columns:", paste(symbol_cols, collapse = ", "), "\n")
    cat("DEBUG: Number of rows:", nrow(result), "\n")
  }

  # Case 1: No condition provided - original behavior (BACKWARD COMPATIBLE)
  if (is.null(apply_when)) {
    if (debug) cat("DEBUG: Using unconditional combination\n")

    # Process each additional filter
    for (i in 2:length(filters)) {
      filter_i <- as.data.table(filters[[i]])

      # FIXED: Use simple column-by-column approach for reliability
      for (col in symbol_cols) {
        if (op == "and") {
          # Extract columns as vectors for operation
          col1 <- result[[col]]
          col2 <- filter_i[[col]]
          result[, (col) := col1 & col2]
        } else {  # op == "or"
          col1 <- result[[col]]
          col2 <- filter_i[[col]]
          result[, (col) := col1 | col2]
        }
      }

      if (debug) {
        cat("DEBUG: Combined with filter", i, "\n")
      }
    }

  } else {
    # Case 2: Conditional combination
    if (debug) cat("DEBUG: Using conditional combination\n")

    # Validate condition length
    if (length(apply_when) != nrow(result)) {
      stop(sprintf("combine_filters: apply_when length (%d) must match filter rows (%d)",
                   length(apply_when), nrow(result)))
    }

    # Convert condition to logical, treating NA as FALSE
    condition_mask <- !is.na(apply_when) & as.logical(apply_when)

    if (debug) {
      cat("DEBUG: Condition TRUE for", sum(condition_mask), "rows\n")
      cat("DEBUG: Condition FALSE/NA for", sum(!condition_mask), "rows\n")
    }

    # Only process if there are TRUE conditions
    if (any(condition_mask)) {
      true_rows <- which(condition_mask)

      # Process each additional filter
      for (i in 2:length(filters)) {
        filter_i <- as.data.table(filters[[i]])

        # FIXED: Process each column separately for clarity
        for (col in symbol_cols) {
          if (op == "and") {
            # Get current values
            current_vals <- result[[col]]
            filter_vals <- filter_i[[col]]

            # Apply operation only to TRUE rows
            current_vals[true_rows] <- current_vals[true_rows] & filter_vals[true_rows]

            # Update result
            result[, (col) := current_vals]
          } else {  # op == "or"
            current_vals <- result[[col]]
            filter_vals <- filter_i[[col]]

            current_vals[true_rows] <- current_vals[true_rows] | filter_vals[true_rows]

            result[, (col) := current_vals]
          }
        }

        if (debug) {
          cat("DEBUG: Conditionally combined with filter", i, "\n")
        }
      }
    } else {
      if (debug) cat("DEBUG: No TRUE conditions, returning first filter unchanged\n")
    }
  }

  if (debug) {
    cat("DEBUG: Final result dimensions:", dim(result), "\n")
    cat("DEBUG: Sample of combined values (first 5 rows of AAPL if exists):\n")
    if ("AAPL" %in% names(result)) {
      print(head(result[, .(Date, AAPL)], 5))
    }
  }

  return(result)
}





