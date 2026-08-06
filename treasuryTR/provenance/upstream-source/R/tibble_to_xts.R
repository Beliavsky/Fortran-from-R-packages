#' Convert a tibble data frame to an xts object
#'
#' @param x a tibble with a column 'date'
#'
#' @return
#' An xts object
#'
#' @export
tibble_to_xts <- function(x) {
  if (!is.data.frame(x)) {
    stop("x must be a (tibble) data frame")
  }
  if (!("date" %in% names(x))) {
    stop("x must have a column `date`")
  }
  if (!(is.Date(x$date))) {
    stop("column `date` must be of date format")
  }

  out <- xts(x = x %>% select(-date),
             order.by = x$date)
  return(out)
}
