#' @title Analytically - Tractable First Passage (AT1P) model
#'
#' @description
#' \code{at1p} calculates the survival probability \eqn{Q(\tau > t)} and default intensity
#' for each maturity according to the structural Analytically - Tractable First Passage model.
#'
#' @param V0 firm value at time \code{t = 0} (it is a constant value).
#' @param H0 value of the safety level at time \code{t = 0}.
#' @param B free positive parameter used for shaping the barrier \code{Ht}.
#' @param sigma a vector of constant stepwise volatility \eqn{\sigma_t}.
#' @param r a vector of constant stepwise risk-free rate.
#' @param t a vector of debt maturity structure (it is a numeric vector).
#'
#' @details
#' In this function the safety level \code{Ht} is calculated using the formula:
#' \deqn{H(t) = \frac{H0}{V0} * E_0[V_t] * \exp^{- B \int_0^t \sigma_u du}}
#' The backbone of the default barrier at \eqn{t} is a proportion, controlled by the parameter
#' \code{H0}, of the expected value of the company assets at \eqn{t}. \code{H0} may depend on the
#' level of the liabilities, on safety covenants, and in general on the characteristics of the capital
#' structure of the company. Also, depending on the value of the parameter \code{B}, it is possible
#' that this backbone is modified by accounting for the volatility of the company assets. For
#' example, if \code{B > 0} corresponds to the interpretation that when volatility increases - which
#' can be independed of credit quality - the barrier is slightly lowered to cut some more slack
#' to the company before going bankrupt. When \code{B = 0} the barrier does not depend on the
#' volatility and the "distance to default" is simply modelled through the barrier parameter \code{H0}.
#'
#' @return
#' \code{at1p} returns an object of class \code{data.frame} containing the firm value, safety level \eqn{H(t)}
#' and the survival probability for each maturity. The last column is the default intensity calculated
#' among each interval \eqn{\Delta t}.
#'
#' @references
#' Damiano Brigo, Massimo Morini, Andrea Pallavicini (2013)
#' Counterparty Credit Risk, Collateral and Funding.
#' With Pricing Cases for All Asset Classes.
#'
#' @examples
#' mod <- at1p(V0 = 1, H0 = 0.7, B = 0.4, sigma = rep(0.1, 10), r = cdsdata$ED.Zero.Curve,
#' t = cdsdata$Maturity)
#' mod
#'
#' plot(cdsdata$Maturity, mod$Ht, type = 'b', xlab = 'Maturity', ylab = 'Safety Level H(t)',
#' main = 'Safety level for different maturities', ylim = c(min(mod$Ht), 1.5),
#' col = 'red')
#' lines(cdsdata$Maturity, mod$Vt, xlab = 'Maturity', ylab = 'V(t)',
#' main = 'Value of the Firm \n at time t', type = 's')
#'
#' plot(cdsdata$Maturity, mod$Survival, type = 'b',
#' main = 'Survival Probability for different Maturity \n (AT1P model)',
#' xlab = 'Maturity', ylab = 'Survival Probability')
#'
#' matplot(cdsdata$Maturity, mod$Default.Intensity, type = 'l', xlab = 'Maturity',
#' ylab = 'Default Intensity')
#'
#'@export
#'
at1p <- function(V0, H0, B, sigma, r, t){

  if(!is.numeric(V0)) stop("V0 must be a number")
  if(length(V0) != 1) stop("V0 must be a number and not a vector")
  if(!is.numeric(H0)) stop("H0 must be a number")
  if(length(H0) != 1) stop("H0 must be a number and not a vector")
  if(!is.numeric(B)) stop("B must be a number greater or equal to zero")
  if(B < 0) stop("B must be a number greater or equal to zero")
  if(length(B) != 1) stop("B must be a number and not a vector")
  if(!is.numeric(sigma)) stop("sigma must be a numeric vector")
  if(!is.numeric(r)) stop("r must be a numeric vector")
  if(!is.numeric(t)) stop("Time t must be a numeric vector")
  if(length(t) != length(sigma)) stop("sigma and t must have the same length")
  if(length(r) != length(t)) stop("r and t must have the same length")

	dt <- diff(c(0, t))

	sigma2u <- cumsum(sigma^2 * dt)

	#Survival probability
	A1 <- stats::pnorm((log(V0 / H0) + 0.5 * (2 * B - 1) * sigma2u) / sqrt(sigma2u))
	A2 <- (H0 / V0)^(2 * B - 1)
	A3 <- stats::pnorm((log(H0 / V0) + 0.5 * (2 * B - 1) * sigma2u) / sqrt(sigma2u))

	Surv <- A1 - (A2 * A3)
	Surv[is.nan(Surv)] <- A1[is.nan(Surv)]

  #Firms value
	Vt <- V0 * exp(r * t)

	#Safety level H(t)
	Ht <- rep(0, length(t))
	for(i in 1 : length(t)){
		Ht[i] <- (H0 / V0) * Vt[i] * exp(- B * sigma2u[i])
	}

	#Default Intensity
	int <- rep(0, length(t))
	int[1] <- -log(Surv[1]) / dt[1]
	for(i in 2 : length(t)){
		int[i] <- (-log(Surv[i]) + log(Surv[i - 1]) ) / dt[i]
    }

	out <- data.frame(t, Vt, Ht, Surv, int)
	names(out) <- c('Maturity', 'Vt', 'Ht', 'Survival', 'Default.Intensity')

	return(out)
}
