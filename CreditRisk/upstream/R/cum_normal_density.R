#' Cumulative Normal Distribution Function
#'
#' This function calculates the cumulative normal distribution function (CDF) 
#' for a given value x using the Hastings approximation method. This approximation 
#' is typically used in finance for the calculation of option pricing probabilities.
#'
#' @param x A numeric value or vector for which the cumulative normal distribution 
#'   is to be calculated.
#'
#' @details
#' The function uses a polynomial approximation as described by E.G. Haug in 
#' "The Complete Guide to Option Pricing Formulas" to estimate the CDF of a normal distribution.
#' The coefficients used in the approximation are specifically chosen to minimize the error
#' in the tail of the distribution, which is critical for financial applications like option pricing.
#'
#' The polynomial approximation is applied to the normal density function:
#' \deqn{N(x) = \frac{1}{\sqrt{2\pi}} e^{-x^2/2}}
#'
#' Then, the cumulative probability is adjusted based on the sign of x:
#' - If x is non-negative, it returns \(1 - t\), where t is the polynomial approximation.
#' - If x is negative, it returns \(t\).
#'
#' The cumulative normal distribution function is important in statistics for hypothesis
#' testing and in finance for the Black-Scholes option pricing formula.
#'
#' @return
#' Returns the cumulative probability under the normal curve from \(-\deqn{\infty}{infinity}\) to x.
#'
#' @references
#' Haug, E.G., The Complete Guide to Option Pricing Formulas.
#' Hastings, C. Approximations for Digital Computers. Princeton Univ. Press, 1955.
#'
#' @examples
#' cum_normal_density(1.96)  
#' cum_normal_density(-1.96) 
#'
#' @export
#' 
cum_normal_density <- function(x) {
  k <- 1 / (1 + 0.2316419 * abs(x))
  a1 <- 0.319381530
  a2 <- -0.356563782
  a3 <- 1.781477937
  a4 <- -1.821255978
  a5 <- 1.330274429
  
  # Calculate the normal density function
  normal_density <- exp(-x^2 / 2) / sqrt(2 * pi)
  
  # Applying the polynomial approximation for the tail
  t <- normal_density * (a1 * k + a2 * k^2 + a3 * k^3 + a4 * k^4 + a5 * k^5)
  
  # Final adjustment for cumulative probability using ifelse for vectorization
  result <- ifelse(x >= 0, 1 - t, t)
  
  # Return the result
  return(result)
}
