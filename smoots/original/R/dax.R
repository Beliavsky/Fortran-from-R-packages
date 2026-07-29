#' German Stock Market Index (DAX) Financial Time Series Data
#'
#' A dataset that contains the daily financial data of the DAX from
#' 1990 to July 2019 (currency in EUR).
#'
#' @format A data frame with 7475 rows and 9 variables:
#' \describe{
#'   \item{Year}{the observation year}
#'   \item{Month}{the observation month}
#'   \item{Day}{the observation day}
#'   \item{Open}{the opening price of the day}
#'   \item{High}{the highest price of the day}
#'   \item{Low}{the lowest price of the day}
#'   \item{Close}{the closing price of the day}
#'   \item{AdjClose}{the adjusted closing price of the day}
#'   \item{Volume}{the traded volume}
#' }
#'@source The data was obtained from Yahoo Finance (accessed: 2019-08-22).
#'
#'\url{https://query1.finance.yahoo.com/v7/finance/download/^GDAXI?period1=631148400&period2=1564524000&interval=1d&events=history&crumb=Iaq1EPZAQRb}
"dax"
