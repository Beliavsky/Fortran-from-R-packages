# globals.R
# Global variables declaration to avoid R CMD check NOTEs
# These are primarily data.table special variables and column names

# Suppress R CMD check notes about undefined global variables
utils::globalVariables(c(
  # data.table special variables
  ".",
  ".SD",
  ".N",
  ":=",
  ".I",
  ".GRP",
  ".BY",

  # data.table .. variables (used in column selection)
  "..col_order",
  "..sector_symbols",
  "..selected_symbols",
  "..symbol_cols",
  "..weight_cols",

  # Common column names used across functions
  "Date",
  "Symbol",
  "Shares",
  "AAPL",  # Used in examples

  # Sector-related columns
  "Industry",
  "SubIndustry",
  "Sector",

  # Breadth analysis columns
  "Breadth_Market",
  "Breadth_Percent",

  # Benchmark columns
  "benchmark_price",
  "benchmark_return",

  # Time-related columns created by data.table
  "year",
  "quarter",
  "year_month",
  "week_end",
  "period",

  # Value columns
  "value",

  # data.table functions that appear as variables
  "first",
  "last",

  # Pipe operator (if using magrittr)
  "%>%",
  # data.table NSE helpers used in ML utilities
  "..syms",
  "score",
  "sharpe",
  "sharpe_med",
  "annret_med",
  #for vol target
  'CASH'

))
utils::globalVariables(c("..sy", "w_lag", "w", "turnover", "Slice", "Ret", "Ret_net"))

# Note: The following are functions that should be imported via NAMESPACE,
# not declared as global variables:
# - RSQLite functions: dbConnect, SQLite, dbDisconnect, dbGetQuery, dbExecute
# - quantmod functions: getSymbols, index, Cl
# - TTR functions: RSI, CCI
# - stats functions: IQR, acf, approx, etc.
# - graphics functions: barplot, hist, legend, text
# - utils functions: write.csv
