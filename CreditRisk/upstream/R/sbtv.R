#' @title Scenario Barrier Time-Varying Volatility AT1P model
#'
#' @description
#' \code{sbtv} calculates the survival probability \eqn{Q(\tau > t)} and default intensity
#' for each maturity according to the structural SBTV model.
#'
#' @param V0 firm value at time \code{t = 0} (it is a constant value).
#' @param H vector of differents safety level at time \code{t = 0}.
#' @param p vector of the probability of different scenario (sum of p must be 1).
#' @param B free positive parameter used for shaping the barrier \code{Ht}.
#' @param sigma a vector of constant stepwise volatility \eqn{\sigma_t}.
#' @param r a vector of constant stepwise risk-free rate.
#' @param t a vector of debt maturity structure (it is a numeric vector).
#'
#' @details
#' \code{sbtv} is an extension of the \code{at1p} model. In this model the parameter \code{H0} used
#' in the \code{at1p} model is replaced by a random variable assuming different values in different
#' scenarios, each scenario with a different probability. The survival probability is calculated
#' as a weighted avarage of the survival probability using the formula:
#' \deqn{SBTV.Surv = \sum_{i = 1}^N p[i] * AT1P.Surv(H[i])}
#' where \code{AT1P.Surv(H[i])} is the survival probability computed according to the AT1P model
#' when \eqn{H_0 = H[i]} and with weights equal to the probabilities of the different scenarios.
#'
#' @return
#' \code{sbtv} returns an object of class \code{data.frame} containing the survival probability
#' for each maturity. The last column is the default intensity calculated
#' among each interval \eqn{\Delta t}.
#'
#' @references
#' Damiano Brigo, Massimo Morini, Andrea Pallavicini (2013)
#' Counterparty Credit Risk, Collateral and Funding.
#' With Pricing Cases for All Asset Classes.
#'
#' @examples
#' mod <- sbtv(V0 = 1, H = c(0.4, 0.8), p = c(0.95, 0.05), B = 0, sigma = rep(0.20, 10),
#'             r = cdsdata$ED.Zero.Curve, t = cdsdata$Maturity)
#' mod
#'
#' plot(cdsdata$Maturity, mod$Survival, type = 'b')
#'
#' @export
#'
sbtv <- function(V0, H, p, B, sigma, r, t){

  if(!is.numeric(V0)) stop("V0 must be a number")
  if(length(V0) != 1) stop("V0 must be a number and not a vector")
  if(!is.numeric(H)) stop("H must be a numeric vector")
  if(!is.numeric(p)) stop("p must be a numeric vector of probability")
  if(sum(p) != 1) stop("the sum of p must be 1 because they are the probability of different scenario")
  if(length(H) != length(p)) stop("H and p must have the same length")
  if(!is.numeric(B)) stop("B must be a number greater or equal to zero")
  if(B < 0) stop("B must be a number greater or equal to zero")
  if(length(B) != 1) stop("B must be a number and not a vector")
  if(!is.numeric(sigma)) stop("sigma must be a numeric vector")
  if(!is.numeric(r)) stop("r must be a numeric vector")
  if(!is.numeric(t)) stop("Time t must be a numeric vector")
  if(length(t) != length(sigma)) stop("sigma and t must have the same length")

	dt <- diff(c(0, t))

	sigma2u <- cumsum(sigma^2 * dt)

	#Survival probability
	Surv <- rep(0, length(t))
	for(i in 1 : length(H)){
		tmp <- p[i] * at1p(V0 = V0, H0 = H[i], B = B, sigma = sigma, r = r, t = t)$Survival
		Surv <- Surv + tmp
	}

	#Default Intensity
	int <- rep(0, length(t))
	int[1] <- -log(Surv[1]) / dt[1]
	for(i in 2 : length(t)){
		int[i] <- (-log(Surv[i]) + log(Surv[i - 1])) / dt[i]
    }

    out <- data.frame(t, Surv, int)
    names(out) <- c('Maturity', 'Survival', 'Default.Intensity')
    return(out)
}
