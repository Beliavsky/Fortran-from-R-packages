#' @title Merton's model
#'
#' @description
#' \code{Merton} calculates the survival probability \eqn{Q(\tau > T)} for
#' each maturity according to the structural Merton's model.
#'
#' @param L debt face value at maturity \code{t = T}; if the value of
#' the firm \eqn{V_T} is below the debt face value to be paid in \eqn{T} the
#' company default has occurred (it is a constant value).
#' @param V0 firm value at time \code{t = 0} (it is a constant value).
#' @param sigma volatility (constant for all t).
#' @param r risk-free rate (constant for all t).
#' @param t a vector of debt maturity structure. The last value of this
#' vector rapresents the debt maturity T.
#'
#' @details
#' In Merton model the default event can occur only at debt maturity T and not before.
#' In this model the debt face value \code{L} represents the constant safety
#' level. In this model the firm value is the sum of the firm equity value \code{St} and
#' ad the firm debt value \code{Dt}. The debt value at time \eqn{t < T} is calculated by the formula:
#' \deqn{D_t = L * \exp(-r  (T - t)) - Put(t, T; V_t, L)}
#' The equity value can be derived as a difference between the firm value and the debt:
#' \deqn{S_t = V_t - D_t = V_t - L * \exp(-r  (T - t)) + Put(t, T; V_t, L) = Call(t, T; V_t, L)}
#' (by the put-call parity) so that in the Merton model the equity can be interpreted as a
#' Call option on the value of the firm.
#'
#' @return
#' \code{Merton} returns an object of class \code{data.frame} with:
#' \itemize{
#'     \item \code{Vt}: expected Firm value at time \eqn{t < T} calculated by the simple formula
#'     \eqn{V_t = V_0 * \exp(r t)}.
#'     \item \code{St}: firm equity value at each \eqn{t < T}. This value can be seen as a call
#'     option on the firm value \code{V_t}.
#'     \item \code{Dt}: firm debt value at each \eqn{t < T}.
#'     \item \code{Survival}: survival probability for each maturity.
#' }
#'
#' @references
#' Damiano Brigo, Massimo Morini, Andrea Pallavicini (2013)
#' Counterparty Credit Risk, Collateral and Funding.
#' With Pricing Cases for All Asset Classes
#'
#' @examples
#' mod <- Merton(L = 10, V0 = 20, sigma = 0.2, r = 0.005,
#'               t = c(0.50, 1.00, 2.00, 3.25, 5.00, 10.00, 15.00, 20.00))
#' mod
#'
#' plot(c(0.50, 1.00, 2.00, 3.25, 5.00, 10.00, 15.00, 20.00), mod$Surv,
#'      main = 'Survival Probability for different Maturity \n (Merton model)',
#'      xlab = 'Maturity', ylab = 'Survival Probability', type = 'b')
#' @export
#'
Merton <- function(L, V0, sigma, r, t){
  
  if(!is.numeric(L)) stop("L must be a number")
  if(length(L) != 1) stop("L must be a number and not a vector")
  if(!is.numeric(V0)) stop("V0 must be a number")
  if(length(V0) != 1) stop("V0 must be a number and not a vector")
  if(!is.numeric(sigma)) stop("sigma must be a number")
  if(length(sigma) != 1) stop("sigma must be a number and not a vector")
  if(!is.numeric(r)) stop("r must be a number")
  if(length(r) != 1) stop("r must be constant for all maturities")
  if(!is.numeric(t)) stop("Time t must be a numeric vector")
  
  d1 <- (log(V0 / L) + (r + (sigma^2 / 2) * t)) / (sigma * sqrt(t))
  d2 <- d1 - sigma * sqrt(t)
  
  #Default Probability
  Q <- stats::pnorm(-d2)
  
  #Survival Probability Q(tau > T)
  Surv <- 1 - Q
  
  #Firm value
  Vt <- V0 * exp(r * t)
  
  #Equity value
  St <- generalized_black_scholes(TypeFlag = "c", S = Vt, X = L, Time = (t[length(t)] - t), r = r, b = r, sigma = sigma)
  
  #Debt value
  Dt <- L * exp(- r * (t[length(t)] - t)) - generalized_black_scholes(TypeFlag = "p", S = Vt, X = L, Time = (t[length(t)] - t), r = r, b = r, sigma = sigma)
  
  out <- data.frame(t, Vt, St, Dt, Surv)
  names(out) <- c('Maturity', 'Vt', 'St', 'Dt', 'Survival')
  
  return(out)
} 
  