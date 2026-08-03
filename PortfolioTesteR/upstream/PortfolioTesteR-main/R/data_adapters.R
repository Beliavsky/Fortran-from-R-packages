## data_adapters.R
# Flexible data adapters for momentum strategies library

# Required libraries for all adapters
#library(data.table)

#' Standardize Data to Library Format
#'
#' @description
#' Internal function that converts various data formats to standard
#' wide format with Date column and symbol columns.
#'
#' @param data Input data in long or wide format
#' @param date_col Name of date column
#' @param symbol_col Name of symbol column (for long format)
#' @param price_col Name of price column (for long format)
#' @param symbol_order Optional ordering for symbols
#'
#' @return Standardized data.table
#' @keywords internal
standardize_data_format <- function(data, date_col = "date", symbol_col = "symbol", price_col = "price", symbol_order = NULL) {
  # Convert to data.table if not already
  setDT(data)

  # Ensure date column is properly formatted
  if (!inherits(data[[date_col]], "Date")) {
    data[, (date_col) := as.Date(get(date_col))]
  }

  # Validate required columns exist
  required_cols <- c(date_col, symbol_col, price_col)
  missing_cols <- required_cols[!required_cols %in% names(data)]
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  # Check which symbols actually have data
  available_symbols <- unique(data[[symbol_col]])
  if (!is.null(symbol_order)) {
    missing_symbols <- symbol_order[!symbol_order %in% available_symbols]
    if (length(missing_symbols) > 0) {
      warning(paste("Symbols not found in data:", paste(missing_symbols, collapse = ", ")))
    }
  }

  # Pivot to wide format: Date, Symbol1, Symbol2, Symbol3, ...
  wide_data <- dcast(data,
                     formula = as.formula(paste(date_col, "~ ", symbol_col)),
                     value.var = price_col)

  # Ensure Date column is named consistently
  setnames(wide_data, old = date_col, new = "Date")

  # Reorder columns to match requested symbol order (if provided)
  if (!is.null(symbol_order)) {
    # Get symbols that actually exist in the data
    existing_symbols <- symbol_order[symbol_order %in% names(wide_data)]
    # Reorder columns: Date first, then symbols in requested order
    col_order <- c("Date", existing_symbols)
    wide_data <- wide_data[, ..col_order]
  }

  # Sort by date
  setorder(wide_data, Date)

  return(wide_data)
}

#' Load Price Data from SQL Database
#'
#' @description
#' Loads stock price data from SQLite database with automatic frequency conversion.
#'
#' @param db_path Path to SQLite database file
#' @param symbols Character vector of stock symbols to load
#' @param start_date Start date (YYYY-MM-DD) or NULL
#' @param end_date End date (YYYY-MM-DD) or NULL
#' @param auto_update Auto-update database before loading (default: TRUE)
#' @param frequency "daily", "weekly", or "monthly" (default: "daily")
#'
#' @return data.table with Date column and one column per symbol
#' @export
#' @examplesIf requireNamespace("RSQLite", quietly = TRUE) && file.exists("sp500.db")
#' \donttest{
#' prices <- sql_adapter(
#'   db_path   = "sp500.db",
#'   symbols   = c("AAPL", "MSFT"),
#'   start_date = "2020-01-01",
#'   end_date   = "2020-12-31",
#'   frequency  = "weekly"
#' )
#' head(prices)
#' }
sql_adapter <- function(db_path, symbols, start_date = NULL, end_date = NULL,
                        auto_update = TRUE, frequency = "daily") {

  # Check if RSQLite is available
  if (!requireNamespace("RSQLite", quietly = TRUE)) {
    stop("RSQLite package required for sql_adapter. Install with: install.packages('RSQLite')")
  }

  #library(RSQLite)

  # Connect to database
  db <- RSQLite::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(RSQLite::dbDisconnect(db))

  # Initialize database if needed
  init_sql_db(db)

  # Auto-update if requested
  if (auto_update) {
    update_symbols_in_db(db, symbols)
  }

  # Build date filter
  date_filter <- ""
  if (!is.null(start_date) || !is.null(end_date)) {
    conditions <- c()
    if (!is.null(start_date)) conditions <- c(conditions, paste("date >= '", start_date, "'", sep=""))
    if (!is.null(end_date)) conditions <- c(conditions, paste("date <= '", end_date, "'", sep=""))
    date_filter <- paste("AND", paste(conditions, collapse = " AND "))
  }

  # Query data
  symbols_str <- paste0("'", paste(symbols, collapse = "','"), "'")
  query <- sprintf("SELECT date, symbol, price FROM stocks WHERE symbol IN (%s) %s ORDER BY date, symbol",
                   symbols_str, date_filter)

  data <- RSQLite::dbGetQuery(db, query)

  if (nrow(data) == 0) {
    stop("No data found for specified symbols and date range")
  }

  # Standardize format with symbol order preserved
  standardized <- standardize_data_format(data, symbol_order = symbols)

  # Convert to weekly if requested
  # Convert to requested frequency
  if (frequency == "weekly") {
    standardized <- convert_to_nweeks(standardized, n = 1)
  } else if (frequency == "monthly") {
    # BUG FIX: Add monthly handling
    standardized <- convert_to_nweeks(standardized, n = 4)
  }

  return(standardized)
}

# Helper function to initialize SQL database
init_sql_db <- function(db) {
  RSQLite::dbExecute(db, "CREATE TABLE IF NOT EXISTS stocks (
    date DATE,
    symbol TEXT,
    price REAL,
    PRIMARY KEY (date, symbol)
  )")

  RSQLite::dbExecute(db, "CREATE INDEX IF NOT EXISTS idx_symbol_date ON stocks(symbol, date)")
}

# Helper function to update symbols in database
update_symbols_in_db <- function(db, symbols) {
  if (!requireNamespace("quantmod", quietly = TRUE)) {
    warning("quantmod package required for auto-update. Skipping update.")
    return()
  }

#  library(quantmod)

  for (symbol in symbols) {
    # Check if symbol needs update
    last_date <- RSQLite::dbGetQuery(db, sprintf("SELECT MAX(date) as last_date FROM stocks WHERE symbol = '%s'", symbol))$last_date

    if (is.na(last_date) || as.Date(last_date) < Sys.Date() - 1) {
      cat("Updating", symbol, "...\n")

      tryCatch({
        from_date <- if(is.na(last_date)) "2014-01-01" else as.Date(last_date)

        stock <- quantmod::getSymbols(symbol, src='yahoo', from=from_date, auto.assign=FALSE, warnings=FALSE)

        if (!is.null(stock)) {
          stock_data <- data.table(
            date = format(as.Date(index(stock)), "%Y-%m-%d"),
            symbol = symbol,
            price = as.numeric(quantmod::Cl(stock))
          )

          # Insert data
          for (i in seq_len(nrow(stock_data))) {
            RSQLite::dbExecute(db, "INSERT OR IGNORE INTO stocks (date, symbol, price) VALUES (?, ?, ?)",
                      params = list(stock_data$date[i], stock_data$symbol[i], stock_data$price[i]))
          }
        }

        Sys.sleep(0.5)  # Rate limiting

      }, error = function(e) {
        warning(paste("Failed to update", symbol, ":", e$message))
      })
    }
  }
}

#' Download Price Data from Yahoo Finance
#'
#' @description
#' Downloads stock price data directly from Yahoo Finance using quantmod.
#' No database required - perfect for quick analysis and experimentation.
#' Get started with real data in under 5 minutes.
#'
#' @param symbols Character vector of stock symbols
#' @param start_date Start date in "YYYY-MM-DD" format
#' @param end_date End date in "YYYY-MM-DD" format
#' @param frequency "daily" or "weekly" (default: "daily")
#'
#' @return Data.table with Date column and one column per symbol
#' @export
#' @examples
#' # Use included sample data
#' data(sample_prices_weekly)
#'
#' # Build a quick momentum strategy with offline data
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
#' selected <- filter_top_n(momentum, n = 2)
#' weights <- weight_equally(selected)
#' result <- run_backtest(sample_prices_weekly, weights, initial_capital = 100000)
#'
#' @examplesIf interactive()
#' # Download tech stocks (requires internet, skipped on CRAN)
#' if (requireNamespace("quantmod", quietly = TRUE)) {
#'   prices <- yahoo_adapter(
#'     symbols = c("AAPL", "MSFT", "GOOGL"),
#'     start_date = "2023-01-01",
#'     end_date = "2023-12-31",
#'     frequency = "weekly"
#'   )
#'   momentum <- calc_momentum(prices, lookback = 12)
#' }
yahoo_adapter <- function(symbols, start_date, end_date, frequency = "daily") {

  if (!requireNamespace("quantmod", quietly = TRUE)) {
    stop("quantmod package required for yahoo_adapter. Install with: install.packages('quantmod')")
  }

 # library(quantmod)

  all_data <- list()

  for (symbol in symbols) {
    cat("Downloading", symbol, "...\n")

    tryCatch({
      stock <- quantmod::getSymbols(symbol, src='yahoo', from=start_date, to=end_date, auto.assign=FALSE, warnings=FALSE)

      if (!is.null(stock)) {
        symbol_data <- data.table(
          date = as.Date(index(stock)),
          symbol = symbol,
          price = as.numeric(quantmod::Cl(stock))
        )
        all_data[[symbol]] <- symbol_data
      }

      Sys.sleep(0.5)  # Rate limiting

    }, error = function(e) {
      warning(paste("Failed to download", symbol, ":", e$message))
    })
  }

  if (length(all_data) == 0) {
    stop("No data downloaded successfully")
  }

  # Combine all data
  combined_data <- rbindlist(all_data)

  # Standardize format with symbol order preserved
  standardized <- standardize_data_format(combined_data, symbol_order = symbols)

  # Convert to weekly if requested
  # Convert to weekly if requested
  if (frequency == "weekly") {
    standardized <- convert_to_nweeks(standardized, n = 1)
  }

  return(standardized)
}

#' Load Price Data from CSV File
#'
#' @description
#' Reads stock price data from CSV files with flexible column naming.
#' Automatically standardizes to library format.
#'
#' @param file_path Path to CSV file
#' @param date_col Name of date column (default: "date")
#' @param symbol_col Name of symbol column (default: "symbol")
#' @param price_col Name of price column (default: "close")
#' @param frequency Target frequency: "daily" or "weekly" (default: "daily")
#' @param symbol_order Optional vector to order symbols
#'
#' @return Data.table with Date column and price columns
#' @export
#' @examples
#' # Create a temporary tidy CSV from included weekly sample data (offline, fast)
#' data("sample_prices_weekly")
#' PW <- as.data.frame(sample_prices_weekly)
#' syms <- setdiff(names(PW), "Date")[1:2]
#'
#' stk  <- stack(PW[1:10, syms])
#' tidy <- data.frame(
#'   Date   = rep(PW$Date[1:10], times = length(syms)),
#'   Symbol = stk$ind,
#'   Price  = stk$values
#' )
#'
#' tmp <- tempfile(fileext = ".csv")
#' write.csv(tidy, tmp, row.names = FALSE)
#' prices <- csv_adapter(tmp)
#' head(prices)
#' unlink(tmp)
csv_adapter <- function(file_path, date_col = "Date", symbol_col = "Symbol", price_col = "Price", frequency = "daily", symbol_order = NULL) {

  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }

  # Read CSV
  data <- fread(file_path)

  # If no symbol_order provided, use the order they appear in the data
  if (is.null(symbol_order)) {
    symbol_order <- unique(data[[symbol_col]])
  }

  # Standardize format with symbol order preserved
  standardized <- standardize_data_format(data, date_col, symbol_col, price_col, symbol_order = symbol_order)

  # Convert to weekly if requested
  if (frequency == "weekly") {
    standardized <- convert_to_nweeks(standardized, n = 1)
  }

  return(standardized)
}

#' Adapter for User-Provided Data
#'
#' @description
#' Simple adapter for when users provide their own data frame.
#' Ensures proper Date formatting and sorting.
#'
#' @param data User-provided data frame
#' @param date_col Name of date column (default: "Date")
#'
#' @return Standardized data.table
#' @export
#' @examples
#' # Use your own data frame
#' data("sample_prices_weekly")
#' my_prices <- manual_adapter(sample_prices_weekly)
manual_adapter <- function(data, date_col = "Date") {
  # For when users provide their own data.frame/data.table
  setDT(data)

  # Ensure Date column
  if (!inherits(data[[date_col]], "Date")) {
    data[, (date_col) := as.Date(get(date_col))]
  }

  setnames(data, old = date_col, new = "Date")
  setorder(data, Date)

  return(data)
}

#' Validate Data Format for Library Functions
#'
#' @description
#' Checks that data meets library requirements: proper Date column,
#' at least one symbol, correct data types. Prints diagnostic info.
#'
#' @param data Data frame to validate
#'
#' @return TRUE if valid, stops with error if not
#' @export
#' @examples
#' data("sample_prices_weekly")
#' # Check if data is properly formatted
#' validate_data_format(sample_prices_weekly)
validate_data_format <- function(data) {
  required_conditions <- list(
    "Must be data.frame or data.table" = is.data.frame(data),
    "Must have Date column" = "Date" %in% names(data),
    "Date column must be Date class" = inherits(data$Date, "Date"),
    "Must have at least one symbol column" = ncol(data) >= 2,
    "Must have at least one row" = nrow(data) >= 1
  )

  failed_conditions <- names(required_conditions)[!unlist(required_conditions)]

  if (length(failed_conditions) > 0) {
    stop("Data validation failed:\n", paste("-", failed_conditions, collapse = "\n"))
  }

  cat("Data validation passed [OK]\n")
  cat("Shape:", nrow(data), "dates x", ncol(data)-1, "symbols\n")
  cat("Date range:", as.character(min(data$Date)), "to", as.character(max(data$Date)), "\n")
  cat("Symbols:", paste(setdiff(names(data), "Date"), collapse = ", "), "\n")

  return(TRUE)
}

# Debug function to check what symbols are available in SQL database
check_sql_symbols <- function(db_path, symbols = NULL, start_date = NULL, end_date = NULL) {
  if (!requireNamespace("RSQLite", quietly = TRUE)) {
    stop("RSQLite package required")
  }

 # library(RSQLite)
  db <- RSQLite::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(RSQLite::dbDisconnect(db))

  # Get all available symbols
  all_symbols <- RSQLite::dbGetQuery(db, "SELECT DISTINCT symbol FROM stocks ORDER BY symbol")$symbol
  cat("All symbols in database:", paste(all_symbols, collapse = ", "), "\n")

  if (!is.null(symbols)) {
    # Check which requested symbols are available
    available <- symbols[symbols %in% all_symbols]
    missing <- symbols[!symbols %in% all_symbols]

    cat("Requested symbols:", paste(symbols, collapse = ", "), "\n")
    cat("Available symbols:", paste(available, collapse = ", "), "\n")
    if (length(missing) > 0) {
      cat("Missing symbols:", paste(missing, collapse = ", "), "\n")
    }

    # Check date ranges for available symbols
    for (symbol in available) {
      date_query <- sprintf("SELECT MIN(date) as min_date, MAX(date) as max_date, COUNT(*) as count FROM stocks WHERE symbol = '%s'", symbol)
      if (!is.null(start_date)) date_query <- paste(date_query, sprintf("AND date >= '%s'", start_date))
      if (!is.null(end_date)) date_query <- paste(date_query, sprintf("AND date <= '%s'", end_date))

      result <- RSQLite::dbGetQuery(db, date_query)
      cat(sprintf("%s: %s to %s (%d records)\n", symbol, result$min_date, result$max_date, result$count))
    }
  }

  return(list(all_symbols = all_symbols, available = available, missing = missing))
}


#' Convert Data to N-Week Frequency
#'
#' @description
#' Resamples daily or weekly data to n-week periods. Handles week-ending
#' calculations and various aggregation methods.
#'
#' @param data Data.table with Date column and price columns
#' @param n Number of weeks to aggregate (default: 1 for weekly)
#' @param method Aggregation method: "last" or "mean" (default: "last")
#'
#' @return Data.table resampled to n-week frequency
#' @export
#' @examples
#' data("sample_prices_daily")
#' # Convert daily to weekly
#' weekly <- convert_to_nweeks(sample_prices_daily, n = 1)
#' # Convert to bi-weekly
#' biweekly <- convert_to_nweeks(sample_prices_daily, n = 2)
convert_to_nweeks <- function(data, n = 1, method = "last") {
  # Convert any frequency to n-week periods
  # n=1 gives weekly, n=2 biweekly, n=4 monthly, etc.

  setDT(data)

  # Ensure Date column is Date class
  if (!inherits(data$Date, "Date")) {
    data[, Date := as.Date(Date)]
  }

  # Get symbol columns
  symbol_cols <- setdiff(names(data), "Date")

  # If n=1 and data is daily, simple weekly conversion
  if (n == 1) {
    # Calculate week ending (Friday) - fixed calculation
    data[, week_end := Date + (5 - as.numeric(format(Date, "%u"))) %% 7]

    if (method == "last") {
      # Include Date in the aggregation to avoid the warning
      weekly_data <- data[, c(list(Date = max(week_end)),
                              lapply(.SD, function(x) tail(na.omit(x), 1))),
                          by = week_end, .SDcols = symbol_cols]
    } else if (method == "mean") {
      weekly_data <- data[, c(list(Date = max(week_end)),
                              lapply(.SD, mean, na.rm = TRUE)),
                          by = week_end, .SDcols = symbol_cols]
    }

    # No need to rename - Date is already set correctly
    weekly_data[, week_end := NULL]  # Remove the grouping column
    setorder(weekly_data, Date)
    return(weekly_data)
  }

  # For n > 1, first convert to weekly if needed
  date_diffs <- diff(data$Date)
  median_diff <- median(as.numeric(date_diffs))

  if (median_diff < 6) {  # Daily data
    # Recursive call with n=1 to get weekly first
    weekly_data <- convert_to_nweeks(data, n = 1, method = method)
  } else {
    weekly_data <- copy(data)  # Already weekly or lower frequency
  }

  # Now resample weekly to n-week periods
  weekly_data[, period := floor((as.numeric(Date - min(Date)) / 7) / n)]

  if (method == "last") {
    nweek_data <- weekly_data[, c(list(Date = max(Date)),
                                  lapply(.SD, function(x) tail(na.omit(x), 1))),
                              by = period, .SDcols = symbol_cols]
  } else if (method == "mean") {
    nweek_data <- weekly_data[, c(list(Date = max(Date)),
                                  lapply(.SD, mean, na.rm = TRUE)),
                              by = period, .SDcols = symbol_cols]
  }

  nweek_data[, period := NULL]
  setorder(nweek_data, Date)

  return(nweek_data)
}


#' Load Adjusted Price Data from SQL Database
#'
#' @description
#' Loads adjusted stock prices (for splits/dividends) from SQLite.
#'
#' @param db_path Path to SQLite database file
#' @param symbols Character vector of stock symbols to load
#' @param start_date Start date (YYYY-MM-DD) or NULL
#' @param end_date End date (YYYY-MM-DD) or NULL
#' @param auto_update Auto-update database (default: FALSE)
#' @param frequency "daily", "weekly", or "monthly" (default: "daily")
#' @param use_adjusted Use adjusted prices if available (default: TRUE)
#'
#' @return data.table with Date column and adjusted prices per symbol
#' @export
#' @examplesIf requireNamespace("RSQLite", quietly = TRUE) && file.exists("sp500.db")
#' \donttest{
#' prices <- sql_adapter_adjusted(
#'   db_path   = "sp500.db",
#'   symbols   = c("AAPL", "MSFT"),
#'   start_date = "2020-01-01",
#'   end_date   = "2020-12-31",
#'   frequency  = "monthly"
#' )
#' head(prices)
#' }
sql_adapter_adjusted <- function(db_path, symbols, start_date = NULL, end_date = NULL,
                                 auto_update = FALSE, frequency = "daily",
                                 use_adjusted = TRUE) {  # New parameter!

  # Check if RSQLite is available
  if (!requireNamespace("RSQLite", quietly = TRUE)) {
    stop("RSQLite package required for sql_adapter. Install with: install.packages('RSQLite')")
  }

  #library(RSQLite)

  # Connect to database
  db <- RSQLite::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(RSQLite::dbDisconnect(db))

  # Auto-update if requested (would need modification for adjusted prices)
  if (auto_update) {
    warning("Auto-update not implemented for adjusted prices. Skipping.")
  }

  # Build date filter
  date_filter <- ""
  if (!is.null(start_date) || !is.null(end_date)) {
    conditions <- c()
    if (!is.null(start_date)) conditions <- c(conditions, paste("date >= '", start_date, "'", sep=""))
    if (!is.null(end_date)) conditions <- c(conditions, paste("date <= '", end_date, "'", sep=""))
    date_filter <- paste("AND", paste(conditions, collapse = " AND "))
  }

  # Choose price column
  price_column <- if(use_adjusted) "price_adjusted" else "price"

  # First check if adjusted column exists
  if(use_adjusted) {
    check_query <- "SELECT COUNT(*) as n FROM pragma_table_info('stocks') WHERE name = 'price_adjusted'"
    has_adjusted <- RSQLite::dbGetQuery(db, check_query)$n > 0

    if(!has_adjusted) {
      warning("price_adjusted column not found. Using regular price column.")
      price_column <- "price"
    }
  }

  # Query data
  symbols_str <- paste0("'", paste(symbols, collapse = "','"), "'")
  query <- sprintf("SELECT date, symbol, %s as price FROM stocks WHERE symbol IN (%s) %s ORDER BY date, symbol",
                   price_column, symbols_str, date_filter)

  data <- RSQLite::dbGetQuery(db, query)

  if (nrow(data) == 0) {
    stop("No data found for specified symbols and date range")
  }

  # Check for NULLs in adjusted prices
  if(use_adjusted && price_column == "price_adjusted") {
    null_count <- sum(is.na(data$price))
    if(null_count > 0) {
      warning(sprintf("Found %d NULL adjusted prices. Consider running update_database_adjusted_prices.R", null_count))
    }
  }

  # Standardize format with symbol order preserved
  standardized <- standardize_data_format(data, symbol_order = symbols)

  # Convert to weekly if requested
  # Convert to requested frequency
  if (frequency == "weekly") {
    standardized <- convert_to_nweeks(standardized, n = 1)
  } else if (frequency == "monthly") {
    # BUG FIX: Add monthly handling
    standardized <- convert_to_nweeks(standardized, n = 4)
  }

  return(standardized)
}



#' Load Mixed Symbols Including VIX
#'
#' @description
#' Handles loading regular stocks and VIX together, with VIX loaded
#' separately without auto-update to avoid issues.
#'
#' @param db_path Path to SQLite database
#' @param symbols Character vector including regular stocks and optionally "VIX"
#' @param start_date Start date for data
#' @param end_date End date for data
#' @param frequency Data frequency (default: "weekly")
#' @param use_adjusted Use adjusted prices (default: TRUE)
#'
#' @return data.table with all symbols properly loaded
#' @export
#' @examplesIf requireNamespace("RSQLite", quietly = TRUE) && file.exists("sp500.db")
#' \donttest{
#' mixed <- load_mixed_symbols(
#'   db_path  = "sp500.db",
#'   symbols  = c("AAPL", "MSFT", "VIX"),
#'   start_date = "2020-01-01",
#'   end_date   = "2020-12-31",
#'   frequency  = "weekly"
#' )
#' head(mixed)
#' }
load_mixed_symbols <- function(db_path, symbols, start_date, end_date,
                               frequency = "weekly", use_adjusted = TRUE) {
  # Separate VIX from other symbols
  vix_included <- "VIX" %in% symbols
  regular_symbols <- symbols[symbols != "VIX"]

  # Load regular symbols with auto-update
  if (length(regular_symbols) > 0) {
    if (use_adjusted) {
      regular_data <- sql_adapter_adjusted(
        db_path = db_path,
        symbols = regular_symbols,
        start_date = start_date,
        end_date = end_date,
        frequency = frequency,
        auto_update = TRUE
      )
    } else {
      regular_data <- sql_adapter(
        db_path = db_path,
        symbols = regular_symbols,
        start_date = start_date,
        end_date = end_date,
        frequency = frequency,
        auto_update = TRUE
      )
    }
  }

  # Load VIX separately without auto-update
  if (vix_included) {
    if (use_adjusted) {
      vix_data <- sql_adapter_adjusted(
        db_path = db_path,
        symbols = "VIX",
        start_date = start_date,
        end_date = end_date,
        frequency = frequency,
        auto_update = FALSE
      )
    } else {
      vix_data <- sql_adapter(
        db_path = db_path,
        symbols = "VIX",
        start_date = start_date,
        end_date = end_date,
        frequency = frequency,
        auto_update = FALSE
      )
    }

    # Combine data
    if (length(regular_symbols) > 0) {
      all_data <- merge(regular_data, vix_data, by = "Date", all = TRUE)
    } else {
      all_data <- vix_data
    }
  } else {
    all_data <- regular_data
  }

  return(all_data)
}

