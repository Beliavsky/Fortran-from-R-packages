#' Sample Weekly Stock Prices
#'
#' Weekly closing prices for 20 stocks from 2017-2019. Data includes
#' major stocks from various sectors and is suitable for demonstrating
#' backtesting and technical analysis functions.
#'
#' @format A data.table with 158 rows and 21 columns:
#' \describe{
#'   \item{Date}{Date object, weekly closing date (typically Friday)}
#'   \item{AAPL}{Apple Inc. adjusted closing price}
#'   \item{AMZN}{Amazon.com Inc. adjusted closing price}
#'   \item{BA}{Boeing Co. adjusted closing price}
#'   \item{BAC}{Bank of America Corp. adjusted closing price}
#'   \item{...}{Additional stock symbols with adjusted closing prices}
#' }
#' @source Yahoo Finance historical data, adjusted for splits and dividends
#' @examples
#' data(sample_prices_weekly)
#' head(sample_prices_weekly)
#' # Calculate momentum
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
#' @usage data(sample_prices_weekly)
#' @docType data
#' @keywords datasets
#' @name sample_prices_weekly
"sample_prices_weekly"

#' Sample Daily Stock Prices
#'
#' Daily closing prices for 20 stocks from 2017-2019. Contains the same
#' symbols as sample_prices_weekly but at daily frequency for more
#' granular analysis and performance calculations.
#'
#' @format A data.table with 754 rows and 21 columns:
#' \describe{
#'   \item{Date}{Date object, trading date}
#'   \item{AAPL}{Apple Inc. adjusted closing price}
#'   \item{AMZN}{Amazon.com Inc. adjusted closing price}
#'   \item{BA}{Boeing Co. adjusted closing price}
#'   \item{BAC}{Bank of America Corp. adjusted closing price}
#'   \item{...}{Additional stock symbols with adjusted closing prices}
#' }
#' @source Yahoo Finance historical data, adjusted for splits and dividends
#' @examples
#' data(sample_prices_daily)
#' head(sample_prices_daily)
#' # Get date range
#' range(sample_prices_daily$Date)
#' @usage data(sample_prices_daily)
#' @docType data
#' @keywords datasets
#' @name sample_prices_daily
"sample_prices_daily"

#' S&P 500 Sector Mappings
#'
#' Sector classifications for the stock symbols in the sample datasets.
#' Note: ETFs (SPY, QQQ, etc.) are not included as they represent
#' indices or sectors themselves rather than individual companies.
#'
#' @format A data.table with 18 rows and 2 columns:
#' \describe{
#'   \item{Symbol}{Character, stock ticker symbol}
#'   \item{Sector}{Character, GICS sector classification}
#' }
#' @source S&P 500 constituent data
#' @examples
#' data(sample_sp500_sectors)
#' head(sample_sp500_sectors)
#' # Count stocks per sector
#' table(sample_sp500_sectors$Sector)
#' @usage data(sample_sp500_sectors)
#' @docType data
#' @keywords datasets
#' @name sample_sp500_sectors
"sample_sp500_sectors"
