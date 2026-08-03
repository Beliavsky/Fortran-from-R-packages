#' Prepare Data for Portfolio Analysis
#'
#' This function prepares the input data by filtering based on a specified date range, removing the date column,
#' and handling missing values. It also generates time-related columns and returns the processed data.
#'
#' @param dat A data frame or matrix where the first column is the date and the remaining columns are the data.
#' @param date_column_index The index of the date column in the input data. Default is 1.
#' @param start_date A character string specifying the start date for filtering the data in 'YYYYMMDD' format.
#' Default is '19990101'.
#' @param end_date A character string specifying the end date for filtering the data in 'YYYYMMDD' format.
#' Default is '20231231'.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{x}{A matrix of the filtered data with missing values handled.}
#'   \item{mon}{A vector of integers representing the number of months from the first date in the data.}
#'   \item{count}{A vector of the number of entries per month.}
#'   \item{Date}{A vector of Date objects representing the filtered dates.}
#' }
#'
#' @examples
#' data <- data.frame(Date = c("19990101", "19990115", "19990201", "19990301", "19990315", "19990401"),
#' Var1 = c(1, 2, -99.99, 4, 5, -99.99),
#' Var2 = c(3, -99.99, 6, 7, 8, 9),
#' Var3 = c(10, 11, 12, 13, -99.99, 15))
#' result <- prepare_data(data, date_column_index = 1, start_date = '19990101', end_date = '19990430')
#' print(result)
#'
#' @export
prepare_data <- function(dat, date_column_index = 1, start_date = '19990101', end_date = '20231231') {
  # Remove the date column and filter data based on the date range
  dat <- dat[,-date_column_index]
  dat <- matrix(unlist(dat), ncol = ncol(dat), byrow = "FALSE")
  data <- dat[dat[,1] >= start_date & dat[,1] <= end_date, ]

  # Prepare date and time-related columns
  Date <- as.Date(as.character(data[,1]), '%Y%m%d')
  mon <- month(Date) + (year(Date) - year(Date)[1]) * 12
  count <- as.numeric(table(mon))

  x <- data[,2:ncol(dat)]

  # Handle missing data
  if (any(x == -99.99)) {
    x <- x[, !apply(x, 2, function(col) any(col == -99.99))]
  }

  return(list(x = x, mon = mon, count = count, Date = Date))
}
