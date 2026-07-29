#' Extract Model Fitted Values
#'
#' Generic function which extracts fitted values from a \code{smoots} class
#' object. Both \code{fitted} and \code{fitted.values} can be called.
#'
#' @param object an object from the \code{smoots} class.
#' @param ... included for consistency with the generic function.
#'
#' @export
#'
#' @return
#' Fitted values extracted from a \code{smoots} class object.
#'
#' @author
#'\itemize{
#'\item Sebastian Letmathe (Scientific Employee) (Department of Economics,
#'Paderborn University), \cr
#'}

fitted.smoots <- function(object, ...){
  object$ye
}


