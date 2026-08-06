#' Load data from FRED
#' @param series The series code as found on https://fred.stlouisfed.org/, see details
#' @param ... Additional parameters handed to quantmod::getSymbols
#' @param format_out xts or tibble
#' @param na_locf replace NA's with last observation
#' @param percent_adjust divide raw data by 100
#'
#'
#' @return The yields data as an xts object
#'
#' @details
#' The function is a wrapper for quantmod::getSymbols().
#'
#' Commonly used constant-maturity yield series are:
#'
#' \itemize{
#' \item{DGS1MO: 1-Month Treasury Constant Maturity Rate}
#' \item{DGS3MO: 3-Month Treasury Constant Maturity Rate}
#' \item{DGS6MO: 6-Month Treasury Constant Maturity Rate}
#' \item{DGS1: 1-Year Treasury Constant Maturity Rate}
#' \item{DGS2: 2-Year Treasury Constant Maturity Rate}
#' \item{DGS3: 3-Year Treasury Constant Maturity Rate}
#' \item{DGS5: 5-Year Treasury Constant Maturity Rate}
#' \item{DGS7: 7-Year Treasury Constant Maturity Rate}
#' \item{DGS10: 10-Year Treasury Constant Maturity Rate}
#' \item{DGS20: 20-Year Treasury Constant Maturity Rate}
#' \item{DGS30: 30-Year Treasury Constant Maturity Rate}
#' }
#'
#' @export
#'
#' @examples
#' # US 1-Month yields
#' yields_us_1m <- get_yields(series = "DGS1MO")
get_yields <- function(series = "DGS10", na_locf = TRUE,
                       percent_adjust = TRUE, format_out = "xts", ...) {
  yields <- tryCatch({
    quantmod::getSymbols(Symbols = series,
                         src = "FRED", auto.assign = FALSE, ...)
    },
    error = function(e) {
      message("An error occurred while accessing the data. Please check your internet connection and try again.")
      return(NULL)
      }
    )

  if (!(format_out %in% c("xts", "tibble"))) {
    stop("format_out has to be one of xts or tibble")
  }

  if (na_locf) {
    yields <- na.locf0(yields)
  }

  if (percent_adjust) {
    yields <- yields/100
  }

  if (format_out == "tibble") {
    yields <- xts_to_tibble(yields)
  }

  return(yields)
}
