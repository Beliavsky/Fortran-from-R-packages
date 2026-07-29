#' Extract Model Residuals
#'
#' Generic function which extracts model residuals from a \code{smoots} class
#' object. Both \code{residuals} and its abbreviated form \code{resid} can be called.
#'
#' @param object an object from the \code{smoots} class.
#' @param ... included for consistency with the generic function.
#'
#' @export
#'
#' @return
#' Residuals extracted from a \code{smoots} class object.
#'
#' @author
#'\itemize{
#'\item Sebastian Letmathe (Scientific Employee) (Department of Economics,
#'Paderborn University), \cr
#'}

residuals.smoots <- function(object, ...){
  object$res
}



