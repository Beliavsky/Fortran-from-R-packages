require(xts)

#' Bitcoin price data 2016-19
#'
#' Time series of Bitcoin closing prices from 31 December 2015 to 31 December 2019 (1044 values).
#' This permits the calculation of 4 calendar years of returns.
#'
#' @docType data
#'
#' @usage data(bitcoin)
#'
#' @format An object of class \code{"xts"}.
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(bitcoin)
#' plot(bitcoin)
#' X <- (diff(log(bitcoin))[-1]) * 100
#' plot(X)
"bitcoin"

#' CPI inflation data 1959-2020
#'
#' Time series of US quarterly CPI (consumer price index) data Q4 1959 to Q4 2020 (245 values) for studying inflation.
#' These data were sourced from the OECD webpage and represent the total `perspective' on inflation, including food and energy.
#' They have been based to have a value of 100 in 2015.
#'
#' @docType data
#'
#' @usage data(cpi)
#'
#' @format An object of class \code{"xts"}.
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(cpi)
#' plot(cpi)
#' X <- (diff(log(cpi))[-1]) * 100
#' plot(X)
"cpi"
