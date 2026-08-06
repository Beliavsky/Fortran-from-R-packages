#' Convert an xts object to a tibble data frame
#'
#' @param x an xts object
#'
#' @note
#' If this function is used outside of the package's functions, i.e. using other data than FRED's, then make sure that the tibble does not contain non-numeric columns. An xts object is a matrix with an index attribute (date) and one cannot mix types in a matrix.
#'
#' @return
#' A tibble data frame with first column 'date'
#'
#' @export
xts_to_tibble <- function(x) {
  if (!is.xts(x)) {
    stop("x has to be an xts object")
  }

  out <- as_tibble(x) %>%
    mutate(date = index(x)) %>%
    select(date, everything())

  return(out)
}
