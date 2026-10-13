#' @title SBTV model calibration to market CDS data
#'
#' @description
#' Compares CDS rates quoted on the market with theoric CDS rates calculeted by the function
#' \code{cds} and looks for the parameters to be used into \code{sbtv}
#' for returning the default intensities corresponding to real market CDS rates performing the
#' minimization of the objective function.
#'
#' @param V0 firm value at time \code{t = 0}.
#' @param p vector of the probability of different scenario (sum of p must be 1).
#' @param cdsrate CDS rates from market.
#' @param r a vector of risk-free rate.
#' @param t a vector of debt maturity structure.
#' @param ... additional parameters used in \code{cds} function.
#'
#' @return
#' This function returns an object of class \code{list} with calculated parameters of
#' \code{sbtv} model and the error occurred in the minimization procedure.
#'
#' @details
#' Inside \code{calibrate.sbtv}, the function \code{objfn} takes the input a
#' vector of parameters and returns the mean error occurred estimating CDS rates with
#' \code{cds} function. The inputs used in \code{cds} are the default intensities calculated by
#' the \code{sbtv} function with the calibrated parameters. In particular the error is
#' calculated as:
#' \deqn{\frac{1}{n}\sum_{i=1}^n (c^{ds}-c^{ds}_{mkt})^2.} This quantity is a function of
#' the default intensities and it is the objective function to be minimized in order to take
#' optimal solutions for intensities.
#'
#' @references
#' Damiano Brigo, Massimo Morini, Andrea Pallavicini (2013)
#' Counterparty Credit Risk, Collateral and Funding.
#' With Pricing Cases for All Asset Classes
#'
#' @examples
#' calibrate.sbtv(V0 = 1, p = c(0.95, 0.05), cdsrate = cdsdata$Par.spread,
#' r = cdsdata$ED.Zero.Curve, t = cdsdata$Maturity)
#'
#'@export
#'
calibrate.sbtv <- function(V0, p, cdsrate, r, t, ...){

  if(!is.numeric(V0)) stop("V0 must be a number")
  if(length(V0) != 1) stop("V0 must be a number and not a vector")
  if(!is.numeric(p)) stop("p must be a numeric vector of probability")
  if(sum(p) != 1) stop("the sum of p must be 1 because they are the probability of different scenario")
  if(!is.numeric(r)) stop("r must be a number")
  if(!is.numeric(cdsrate)) stop("cdsrate must be a numeric vector")
  if(!is.numeric(t)) stop("Time t must be a numeric vector")
  if(length(cdsrate) != length(t)) stop("cdsrate and t must have the same length")

	objfn <- function(par, p, cdsrate, r, t){
		B <- par[1]
		sigma <- par[2 : (1 + length(t))]
		H <- par[(2 + length(t)) : length(par)]

		intensity <- sbtv(V0 = V0, H = H, p = p, B = B, sigma = sigma, r = r, t = t)$Default.Intensity
		rate <- cds(t = t, int = intensity, r = r)$Rate

		err <- mean((rate - cdsrate)^2)
		return(err)
	}

	H <- rep(0.4 * V0, length(p))
	for (i in 2 : length(p)){
		H[i] <- H[i - 1] * (1 + 0.5)
	}

	par <- c(0.02, rep(0.2, length(t)), H)
	out <- stats::optim(par = par, fn = objfn, p = p, cdsrate = cdsrate, t = t, r = r, method = 'L-BFGS-B', lower = 0 * t, control = list(maxit = 5000))

	out <- list(V0, out$par[(2 + length(t)) : length(par)], p, out$par[1], out$par[2 : (1 + length(t))], out$value)
	names(out) <- c('V0', 'H', 'p', 'B', 'Sigma', 'Error')
	return(out)
}
