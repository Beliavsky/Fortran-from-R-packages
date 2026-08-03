#' Laplace distribution
#' 
#' Implements some aspects of Laplace distribution (based on stats package) for 
#' stable random walk simulation.
#'
#' @param n numeric, number of observations.
#' @param x numeric, vector of responses.
#' @param b numeric, the scale parameter, where the variance is 2*b^2.
#' 
#' @return numeric, standard convention is followed: 
#'         d* returns the density, 
#'         p* returns the distribution function, 
#'         q* returns the quantile function, and 
#'         r* generates random deviates.
#'
#' @keywords Laplace
#'
#' @author Stephen H-T. Lihn
#'
#' @importFrom stats dexp
#' @importFrom stats rexp
#' @importFrom stats runif
#'
#' @export rlaplace0
#' @export dlaplace0
#' 
### <======================================================================>
rlaplace0 <- function(n, b=1) stats::rexp(n,1/b)*sign(stats::runif(n)-0.5)
### <---------------------------------------------------------------------->
#' @rdname rlaplace0
dlaplace0 <- function(x, b=1) stats::dexp(abs(x),1/b)/2
### <---------------------------------------------------------------------->
