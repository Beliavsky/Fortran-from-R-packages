#' @title Black and Cox's model
#'
#' @description
#' \code{BlackCox} calculates the survival probability \eqn{Q(\tau > t)} and default intensity
#' for each maturity according to the structural Black and Cox's model.
#'
#' @param L debt face value at maturity \code{t = T} (it is a constant value).
#' @param K positive parameter needed to calculate the safety level.
#' @param V0 firm value at time \code{t = 0} (it is a constant value).
#' @param sigma volatility (constant for all t).
#' @param r risk-free rate (constant for all t).
#' @param gamma interest rate used to discount the safety level \code{Ht} (it is a constant value).
#' @param t a vector of debt maturity structure (it is a numeric vector).
#'
#' @details
#' In Merton's model the default event can occurr only at debt maturity \eqn{T} while
#' in Black and Cox's model the default event can occurr even before.
#' In this model the safety level is given by the output \code{Ht}. Hitting this barrier is
#' considered as an erlier default. Assuming a debt face value of \code{L} at the final
#' maturity that coincides with the safety level in \eqn{t = T}, the safety level in \eqn{t\le T} is the
#' \code{K}, with \eqn{K\le L}, value discounted at back at time \eqn{t} using the interest rate
#' \code{gamma}, obtaining: \deqn{H(t | t\le T) = K * \exp^{- \gamma * (T- t)}}
#' The output parameter \code{Default.Intensity} represents the default intensity of
#' \eqn{\Delta t}. The firm's value \code{Vt} is calculated as in the \code{Merton} function.
#'
#' @return
#' This function returns an object of class \code{data.frame} containing firm value, safety level \eqn{H(t)}
#' and the survival probability for each maturity. The last column is the default intensity calculated
#' among each interval \eqn{\Delta t}.
#'
#' @references
#' David Lando  (2004) Credit risk modeling.
#'
#' Damiano Brigo, Massimo Morini, Andrea Pallavicini (2013)
#' Counterparty Credit Risk, Collateral and Funding.
#' With Pricing Cases for All Asset Classes.
#'
#' @examples
#' mod <- BlackCox(L = 0.55, K = 0.40, V0 = 1, sigma = 0.3, r = 0.05, gamma = 0.04,
#' t = c(0.50, 1.00, 2.00, 5.00, 7.00, 10.00, 20.00, 30.00))
#' mod
#'
#' plot(c(0.50, 1.00, 2.00, 5.00, 7.00, 10.00, 20.00, 30.00), mod$Ht, type = 'b',
#'      xlab = 'Maturity', ylab = 'Safety Level H(t)', main = 'Safety level for different
#'      maturities', ylim = c(min(mod$Ht), 1.5), col = 'red')
#' abline(h = 0.55, col = 'red')
#' lines(c(0.50, 1.00, 2.00, 5.00, 7.00, 10.00, 20.00, 30.00), mod$Vt, xlab = 'Maturity',
#'       ylab = 'V(t)', main = 'Value of the Firm \n at time t', type = 's')
#'
#' plot(c(0.50, 1.00, 2.00, 5.00, 7.00, 10.00, 20.00, 30.00), mod$Survival, type = 'b',
#'      main = 'Survival Probability for different Maturity \n (Black & Cox model)',
#'      xlab = 'Maturity', ylab = 'Survival Probability')
#'
#' matplot(c(0.50, 1.00, 2.00, 5.00, 7.00, 10.00, 20.00, 30.00), mod$Default.Intensity,
#'         type = 'l', xlab = 'Maturity', ylab = 'Default Intensity')
#'
#' @export
#'
BlackCox <- function(L, K = L, V0, sigma, r, gamma, t){

  if(!is.numeric(L)) stop("L must be a number")
  if(length(L) != 1) stop("L must be a number and not a vector")
  if(!is.numeric(K)) stop("K must be a number")
  if(length(K) != 1) stop("K must be a number and not a vector")
  if(K > L) stop("K must be less or equal to L")
  if(!is.numeric(V0)) stop("V0 must be a number")
  if(length(V0) != 1) stop("V0 must be a number and not a vector")
  if(!is.numeric(sigma)) stop("sigma must be a number")
  if(length(sigma) != 1) stop("sigma must be a number and not a vector")
  if(!is.numeric(r)) stop("r must be a number")
  if(length(r) != 1) stop("r must be constant for all maturities")
  if(!is.numeric(gamma)) stop("gamma must be a number")
  if(length(gamma) != 1) stop("gamma must be constant for all maturities")
  if(!is.numeric(t)) stop("Time t must be a numeric vector")

  dt <- diff(c(0, t))
	v <- r - gamma - (0.5 * sigma^2)
	a <- v / sigma^2

	#Firms value
	Vt <- V0 * exp(r * t)

	#Safety level H(t)
	Ht <- rep(L, length(t))
	for(i in (length(t) - 1) : 1){
		Ht[i] <- K * exp(- gamma * (t[length(t)] - t[i]))
	}
	H0 <- K * exp(- gamma * t[length(t)])

	#Survival Probability
	A1 <- stats::pnorm((log(V0 / H0) + (v * t)) / (sigma * sqrt(t)))
  A2 <- (H0 / V0)^(2 * a)
  A3 <- stats::pnorm((log(H0 / V0) + (v * t)) / (sigma * sqrt(t)))

	Surv <- A1 - (A2 * A3)
  Surv[is.nan(Surv)] <- A1[is.nan(Surv)]

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
