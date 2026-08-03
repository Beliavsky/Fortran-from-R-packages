#' Package to compute the log likelihood for an inverse gamma stochastic volatility model and to draw from the exact posterior of the inverse volatilities.
#'
#' @description This package computes the log likelihood for an inverse gamma stochastic volatility model using a closed form expression of the likelihood. 
#' 
#' @docType package
#' @name extra
#' @importFrom Rcpp evalCpp
#' @useDynLib invgamstochvol, .registration=TRUE
NULL
#> NULL