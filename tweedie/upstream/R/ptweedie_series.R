#' @title Series Evaluation for the Tweedie Distribution Function
#' @name ptweedie_series
#' @description
#' Evaluates the distribution function (\acronym{df}) for Tweedie distributions 
#' with \eqn{1 < p < 2}{1 < p < 2}
#' using an infinite series, for given values of the dependent variable \code{y}, 
#' the mean \code{mu}, dispersion \code{phi}, and power parameter \code{power}.
#' \emph{Not usually called by general users}, but can be in the case of evaluation problems.
#'
#' @usage ptweedie_series(q, power, mu, phi, verbose = FALSE, details = FALSE)
#' 
#' @param q vector of quantiles.
#' @param power the power parameter \eqn{p}{power}.
#' @param mu the mean parameter \eqn{\mu}{mu}.
#' @param phi the dispersion parameter \eqn{\phi}{phi}.
#' @param verbose logical; if \code{TRUE}, displays some internal computation details. The default is \code{FALSE}.
#' @param details logical; if \code{TRUE}, returns the value of the distribution function and some details.
#' 
#' @return A numeric vector of densities.
#' 
#' @note
#' The 'exact' values for the inverse Gaussian distribution are not really exact, but evaluated using inverse normal distributions,
#' for which very good numerical approximation are available in R.
#' 
#' @references
#' Dunn, Peter K and Smyth, Gordon K (2005).
#' Series evaluation of Tweedie exponential dispersion model densities
#' \emph{Statistics and Computing},
#' \bold{15}(4). 267--280.
#' \doi{10.1007/s11222-005-4070-y}
#' 
#' @examples
#' # Plot a Tweedie distribution function
#' y <- seq(0.01, 4, length = 50)
#' Fy <- ptweedie_series(y, power = 1.1, mu = 1, phi = 1)
#' plot(y, Fy, type = "l", lwd = 2, ylab = "Distribution function")
#' 
#' @importFrom stats dpois 
#'
#' @keywords distribution
#'
#' @export
ptweedie_series <- function(q, power, mu, phi, verbose = FALSE, details = FALSE) {
  ### NOTE: No notation checks

  # SET UP
  lambda <- mu ^ (2 - power) / ( phi * (2 - power) )
  tau    <- phi * (power - 1) * mu ^ ( power - 1 )
  alpha  <- (2 - power) / (1 - power)
  drop <- 39

  # FIND THE LIMITS ON N, the summation index
  # The *lower* limit on N
  lambda <- max(lambda )
  logfmax <-  -log(lambda)/2
  estlogf <- logfmax
  N <- max( lambda )
  
  while ( ( estlogf > (logfmax - drop) ) & ( N > 1 ) ) {
    N <- max(1, N - 2)
    estlogf <- -lambda + N * ( log(lambda) - log(N) + 1 ) - log(N)/2
  }
  lo.N <- max(1, floor(N) )
  
  
  # The *upper* limit on N
  lambda <- min( lambda )
  logfmax <-  -log(lambda) / 2
  estlogf <- logfmax
  N <- max( lambda )
  
  while ( estlogf > (logfmax - drop) ) {
    N <- N + 1
    estlogf <- -lambda + N * ( log(lambda) - log(N) + 1 ) - log(N)/2
  }
  hi.N <- max( ceiling(N) )
  if (verbose) cat("Summing over", lo.N, "to", hi.N, "\n")
  
  # Add a safety check:
  hi.N <- min(hi.N, 1e6)
  if (hi.N < lo.N) hi.N <- lo.N
  
  # EVALUATE between limits of N
  cdf <- array( dim = length(q), 0 )
  
  lambda <- mu ^ (2 - power) / ( phi * (2 - power) )
  tau    <- phi * (power - 1) * mu ^ ( power - 1 )
  alpha  <- (2 - power) / (1 - power)
  
  
  N_vec <- lo.N : hi.N
  
  # The Poisson weights
  pois_den <- dpois(N_vec, lambda)
  
  # The incomplete Gamma values
  # We want a matrix where rows = q and columns = N
  incgamma_matrix <- outer(2 * q / tau, -2 * alpha * N_vec, stats::pchisq)
  
  # Multiply each column (N) by its Poisson weight and sum across rows
  # %*% is a matrix multiplication that does the weighting and summing in one go
  cdf <- as.vector(incgamma_matrix %*% pois_den)
  
  # for (N in (lo.N : hi.N)) {
  #   # Poisson density
  #   pois.den <- dpois( N, lambda)
  #   
  #   # Incomplete gamma
  #   incgamma.den <- stats::pchisq(2 * q / tau, 
  #                          -2 * alpha * N )
  #   
  #   # What we want
  #   cdf <- cdf + pois.den * incgamma.den
  #   
  # }
  
  cdf <- cdf + exp( -lambda )
  its <- hi.N - lo.N + 1
  
  if (details) {
    return( list( cdf = cdf,
                  iterations = its) )
  } else {
    return(cdf)
  }
  
}



#' @rdname ptweedie_series
#' @export
ptweedie.series <- function(q, power, mu, phi, verbose = FALSE, details = FALSE){ 
  lifecycle::deprecate_warn(when = "3.0.5", 
                            what = "ptweedie.series()", 
                            with = "ptweedie_series()")
  ptweedie_series(q, power, mu, phi, verbose = FALSE, details = FALSE)
}

