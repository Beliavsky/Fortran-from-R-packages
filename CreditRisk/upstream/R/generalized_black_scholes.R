#' Generalized Black-Scholes Option Pricing Model
#'
#' This function calculates the price of a European call or put option using the
#' generalized Black-Scholes formula, which extends the standard model to
#' incorporate a continuous dividend yield.
#'
#' @param TypeFlag A character vector indicating the type of option to be priced,
#' either "c" for call options or "p" for put options.
#' @param S Current stock price (scalar).
#' @param X Strike price of the option (scalar).
#' @param Time Time to expiration of the option (in years).
#' @param r Risk-free interest rate (annualized).
#' @param b Cost of carry rate, b = r - q for a dividend yield q.
#' @param sigma Volatility of the underlying asset (annualized).
#'
#' @details
#' The generalized Black-Scholes formula considers both the risk-free rate and
#' a cost of carry, making it suitable for a wider range of financial instruments,
#' including commodities and currencies with continuous yields.
#'
#' The pricing formula for call and put options is determined by:
#' \deqn{C = S e^{(b-r)T} N(d_1) - X e^{-rT} N(d_2)}
#' \deqn{P = X e^{-rT} N(-d_2) - S e^{(b-r)T} N(-d_1)}
#' where:
#' \deqn{d_1 = \frac{\log(S / X) + (b + \sigma^2 / 2) T}{\sigma \sqrt{T}}}
#' \deqn{d_2 = d_1 - \sigma \sqrt{T}}
#' and \eqn{(N(\cdot))} is the cumulative normal distribution function, estimated
#' by the `cum_normal_density` function.
#'
#' @return
#' Returns the price of the specified option (call or put).
#'
#' @references
#' Haug, E.G., The Complete Guide to Option Pricing Formulas.
#'
#' @examples
#' # Calculate the price of a call option
#' generalized_black_scholes("c", S = 100, X = 100, Time = 1, r = 0.05, b = 0.05, sigma = 0.2)
#' # Calculate the price of a put option
#' generalized_black_scholes("p", S = 100, X = 100, Time = 1, r = 0.05, b = 0.05, sigma = 0.2)
#'
#' @export
#' 
generalized_black_scholes = function(TypeFlag = c("c", "p"), S, X, Time, r, b, sigma){
  
  # Description:
  #   Calculate the Generalized Black-Scholes option
  #   price either for a call or a put option.
  
  # References:
  #   Haug E.G., The Complete Guide to Option Pricing Formulas
  
  # FUNCTION:
  
  # Compute:
  TypeFlag = TypeFlag[1]
  d1 = ( log(S/X) + (b+sigma*sigma/2)*Time ) / (sigma*sqrt(Time))
  d2 = d1 - sigma*sqrt(Time)
  
  # price
  if (TypeFlag == "c")
    result = S*exp((b-r)*Time)*cum_normal_density(d1) - X*exp(-r*Time)*cum_normal_density(d2)
  if (TypeFlag == "p")
    result = X*exp(-r*Time)*cum_normal_density(-d2) - S*exp((b-r)*Time)*cum_normal_density(-d1)
  
  
  # Return Value:
  return(result)
}
