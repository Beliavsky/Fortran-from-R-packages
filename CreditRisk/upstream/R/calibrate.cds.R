#' @title Calibrate the default intensities to market CDS data
#'
#' @description
#' Compares CDS rates quoted on market with theoric CDS rates and looks for default intensities
#' that correspond to real market CDS rates trough a minimization problem of an objective
#' function.
#'
#' @param r interest rates.
#' @param t premiums timetable.
#' @param Tj CDS maturities.
#' @param cdsrate CDS rates from market.
#' @param ... additional parameters used in \code{cds} function.
#'
#' @return
#' returns an object of class \code{list} with calculated intensities and the
#' error occurred in the minimization procedure.
#'
#' @details
#' Inside \code{calibrate.cds}, the function \code{err.cds} takes the input a
#' vector of intensities and return the mean error occurred estimating CDS rates with
#' \code{cds}. In particular such error is calculated as:
#' \deqn{\frac{1}{n}\sum_{i=1}^n (c^{ds}-c^{ds}_{mkt})^2.} This quantity is a function of default intensities and is the our objective
#'   function to be minimized in order to take optimal solutions for intensities.
#'
#' @references
#' David Lando (2004) Credit risk modeling
#'
#' Damiano Brigo, Massimo Morini, Andrea Pallavicini (2013)
#' Counterparty Credit Risk, Collateral and Funding.
#' With Pricing Cases for All Asset Classes
#'
#' @examples
#' calibrate.cds( r = cdsdata$ED.Zero.Curve, t = seq(.5, 30, by = 0.5),
#'                Tj = c(1, 2, 3, 4, 5, 7, 10, 20, 30), cdsrate = cdsdata$Par.spread, RR = 0.4)
#'
#' @export
#'
calibrate.cds <- function( r, t, Tj, cdsrate, ...){

	if( !is.numeric(r) ) stop( "interest rate r must be a numeric vector" )
	if( !is.numeric(t) ) stop( "time t must be a numeric vector" )
	if( !is.numeric(cdsrate) ) stop( "cdsrate must be a numeric vector" )

	err.cds <- function(int){
		new.int <- stats::approx( x = Tj, y = int, xout = t, method = "linear", rule = 2)$y
		cds1 <- cds(t = t, r = r, int = new.int, ...)$Rate
		cds1 <- stats::approx( x = t, y = cds1, xout = Tj, rule = 2)$y
		err <- mean((cds1 - cdsrate)^2)
	}

    int = rep(0.1, length(Tj));

	out <- stats::optim( par = int, fn = err.cds,  lower = 0*t, control = list(maxit = 5000), method = "L-BFGS-B")
	out <- list(out$par,out$value)
	names(out) <- c("min_int", "error")

	return(out)
}

