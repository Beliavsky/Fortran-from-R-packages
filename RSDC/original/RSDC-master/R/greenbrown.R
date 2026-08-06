#' Green vs Brown portfolio dataset
#'
#' Daily returns for a green and a brown portfolios constructed following
#' the equal-weighted 10-90 percentile approach.
#'
#' @format A data frame with 2266 rows and three columns:
#' \describe{
#'   \item{DATE}{Dates ranging from 2014-01-02 to 2022-12-30.}
#'   \item{return_green}{Numeric returns for the green portfolio.}
#'   \item{return_brown}{Numeric returns for the brown portfolio.}
#' }
#'
#' @source Built from the source workbook in \code{data-raw/green-brown-ptf.xlsx}
#'   (not shipped on CRAN); see \code{data-raw/greenbrown.R}.
#'
#' @examples
#' data("greenbrown")
#' str(greenbrown)
#' head(greenbrown)
#'
#' @keywords datasets
#' @docType data
#' @usage data(greenbrown)
"greenbrown"
