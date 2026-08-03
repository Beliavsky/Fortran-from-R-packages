#' German Stock Market Index (DAX) Financial Time Series Data
#'
#' A dataset that contains the daily financial data of the DAX from
#' 2000 to December 2023 (currency in EUR).
#'
#' @format A data frame with 6094 rows and 11 variables:
#' \describe{
#'   \item{ticker}{id of the stock}
#'   \item{ref_date}{date in format YY-MM-DD}
#'   \item{price_open}{opening price (daily)}
#'   \item{price_high}{highest price (daily)}
#'   \item{price_low}{lowest price (daily)}
#'   \item{price_close}{closing price (daily)}
#'   \item{volume}{trading volume}
#'   \item{price_adjusted}{adjusted closing price (daily)}
#'   \item{ret_adjusted_prices}{returns obtained from the adj. closing prices}
#'   \item{ret_closing_prices}{returns obtained from the closing prices}
#'   \item{cumret_adjusted_prices}{accumulated arithmetic/log return for the period}
#' }
#' @source The data was obtained from Yahoo Finance.
"DAX"
