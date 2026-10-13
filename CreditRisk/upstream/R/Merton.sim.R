#' @title Firm value in Merton's model
#'
#' @description
#' With this function we simulate \code{n} trajectories of firm value based on
#' Merton's model.
#'
#' @param V0 firm value at time \code{t = 0}.
#' @param r risk-free interest rate (constant for all t).
#' @param sigma volatility (constant for all t).
#' @param t a vector of debt maturity structure.
#' @param n number of trajectories to be generated.
#' @param seed starting seed, default seed is setted randomly.
#'
#' @details
#' The trajectories are calculated according to the equation:
#' \deqn{V_T = V_0 \exp{\int_0^T dln V_t}}
#' Where we express \code{dln V_t} using Ito's lemma to derive the differential
#' of the logarithm of the firm value as:
#' \deqn{dln V_t =(\mu - \frac{\sigma^2}{2})dt + \sigma dW_t}
#'
#' @return
#' This function returns a matrix containing the simulated firm values.
#'
#' @references
#' Gergely Daròczi, Michael Puhle, Edina Berlinger, Péter Csòka, Dàniel Havran
#' Màrton Michaletzky, Zsolt Tulasay, Kata Vàradi, Agnes Vidovics-Dancs (2013)
#' Introduction to R for Quantitative Finance.
#'
#' @examples
#' V <- Merton.sim(V0 = 20, r = 0.05, sigma = 0.2, t = seq(0, 30, by = 0.5), n = 5)
#' matplot(x = seq(0, 30, by = 0.5), y = V, type = 's', lty = 1, xlab = 'Time',
#' ylab = 'Firm value trajectories', main = "Trajectories of the firm values in the Merton's model")
#'
#' @export
#'
Merton.sim <- function(V0, r, sigma, t, n, seed = as.numeric(Sys.time())){

  if(!is.numeric(V0)) stop("V0 must be a number")
  if(length(V0) != 1) stop("V0 must be a number and not a vector")
  if(!is.numeric(r)) stop("r must be a number")
  if(length(r) != 1) stop("r must be constant for all maturities")
  if(!is.numeric(t)) stop("Time t must be a numeric vector")
  if(!is.numeric(n)) stop("n must be a number")
  if(!is.numeric(seed)) stop("seed must be a number")
  if(length(seed) != 1) stop("seed must be a number and not a vector")

	set.seed(seed)

	#Number of time periods
	M <- length(t)

	#length of the period
	dt <- diff(c(0, t))

	#Simulated Values
	val <-  stats::rnorm(n * M, mean = (r - sigma^2 / 2) * dt, sd = sigma * dt^.5)
	dlnV <- matrix(val, M, n)

	#variation of the firm value in time
	V <- V0 * exp(apply(dlnV, 2, cumsum))

	return(V)
}
