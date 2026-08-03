#' The numericOrNull class
#'
#' The S4 class union of numeric and NULL, primarily used for detla
#'
#' @name numericOrNull-class
#'
#'
#' @exportClass numericOrNull
setClassUnion("numericOrNull", c("numeric", "NULL"))

# end

