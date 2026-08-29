#' @title Fourier Inversion Evaluation for the Tweedie Distribution Function
#' @name ptweedie_inversion
#' @description
#' Evaluates the distribution function (\acronym{df}) for Tweedie distributions using Fourier inversion, 
#' for given values of the dependent variable \code{y}, 
#' the mean \code{mu}, dispersion \code{phi}, and power parameter \code{power}.
#' \emph{Not usually called by general users}, but can be in the case of evaluation problems.
#'
#' @usage ptweedie_inversion(q, mu, phi, power, verbose = FALSE, details = FALSE, IGexact = TRUE)
#'
#' @param q vector of quantiles.
#' @param power the power parameter \eqn{p}{power}.
#' @param mu the mean parameter.
#' @param phi the dispersion parameter.
#' @param verbose logical; if \code{TRUE}, displays some internal computation details. The default is \code{FALSE}.
#' @param details logical; if \code{TRUE}, returns the value of the distribution and some information about the integration. The default is \code{FALSE}.
#' @param IGexact logical; if \code{TRUE} (the default), evaluate the inverse Gaussian distribution using the 'exact' values, otherwise uses inversion.
#' 
#' @return If \code{details = FALSE}, a numeric vector of the distribution function values; if \code{details = TRUE}, a list containing \code{CDF} (a vector of the values of the distribution function), \code{regions} (a vector of the number of integration regions used), and \code{exitstatus} (a vector, where a \code{1} for any value means a computational problem or target relative accuracy not reached, for the corresponding observation).
#' 
#' For special cases of \eqn{p} (i.e., \eqn{p = 0, 1, 2, 3}), where no inversion is needed, \code{regions} is set to \code{NA} for all values of \code{q}.
#' For special cases of \code{q} for other values of \eqn{p} (i.e., \eqn{P(Y = 0)}), \code{regions} is set to \code{NA}.
#'
#' @note
#' The 'exact' values for the inverse Gaussian distribution are not really exact, but evaluated using inverse normal distributions,
#' for which very good numerical approximation are available in R.
#' 
#' @references
#' Dunn, P. K. and Smyth, G. K. (2008).
#' Evaluation of Tweedie exponential dispersion model densities by Fourier inversion.
#' \emph{Statistics and Computing}, 
#' \bold{18}, 73--86.
#' \doi{10.1007/s11222-007-9039-6}
#'
#' @examples
#' # Plot a Tweedie distribution function
#' y <- seq(0.01, 4, length = 50)
#' Fy <- ptweedie_inversion(y, mu = 1, phi = 1, power = 1.1)
#' plot(y, Fy, type = "l", lwd = 2, ylab = "Distribution function")
#' 
#' @keywords distribution
#' 
#' @export
ptweedie_inversion <- function(q, mu, phi, power, verbose = FALSE, details = FALSE, IGexact = TRUE ){ 
  ### NOTE: No notation checks
  
  # Check
  if (length(q) == 0L) {
    return(numeric(0))
  }
  
  # CHECK THE INPUTS ARE OK AND OF CORRECT LENGTHS
  if (verbose) cat("- Checking, resizing inputs\n")
  out <- check_inputs(q, mu, phi, power)
  mu <- out$mu
  phi <- out$phi

  # cdf    is the whole vector; the same length as  y.
  # All is resolved in the end.
  cdf <- numeric(length = length(q) )
  regions <- rep(NA, length(q)) 
  
  # IDENTIFY SPECIAL CASES
  special_y_cases <- rep(FALSE, 
                         length(q) )
  if (verbose) cat("- Checking for special cases\n")
  out <- special_cases(q, mu, phi, power,
                       IGexact = IGexact,
                       type = "CDF")
  
  special_p_cases <- out$special_p_cases
  special_y_cases <- out$special_y_cases
  
  if (verbose & special_p_cases) cat("  - Special case for p used\n")
  if ( any(special_y_cases) ) {
    special_y_cases <- out$special_y_cases  
    if (verbose) cat("  - Special cases for first input found\n")
    cdf <- out$f # This is the final vector of results to return, filled with the special-case info.
    # NOTE: regions filled with zeros by default, so regions = 0 in these cases
  }
  
  if ( special_p_cases ) {
    cdf <- out$f
  } else {
    # NOT special p case; ONLY special y cases 
    
    # Now use FORTRAN on the remaining values:
    N_nonSpecial <- length(q) - sum(out$special_y_cases) 
      
    
    
    ### BEGIN SET UP
    pSmall  <- ifelse( (power > 1) & (power < 2),
                       TRUE, 
                       FALSE )
  
    ### END SET UP
  
    if (N_nonSpecial > 0 ) {
      tmp <- .C( "twcomputation",
                 N           = as.integer(N_nonSpecial),              # number of observations
                 power       = as.double(power),                      # p
                 phi         = as.double(phi[!special_y_cases]),      # phi
                 y           = as.double(q[!special_y_cases]),        # y
                 mu          = as.double(mu[!special_y_cases]),       # mu
                 verbose     = as.integer(verbose),                   # verbosity
                 pdf         = as.integer(0),                         # 0: FALSE, as this is the CDF not PDF
                 # THE OUTPUTS:
                 funvalue    = numeric(N_nonSpecial),                 # funvalue
                 exitstatus  = integer(N_nonSpecial),                 # exitstatus
                 relerr      = numeric(N_nonSpecial),                 # relerr
                 its         = integer(N_nonSpecial),                 # its
                 PACKAGE     = "tweedie")
      cdf[!special_y_cases] <- tmp$funvalue
      regions[!special_y_cases] <- tmp$its
    }
  }
  
  if (details) {
    return( list( cdf = cdf,
                  regions = regions,
                  exitstatus = tmp$exitstatus))
  } else {
    return(cdf)
  }
}

#' @rdname ptweedie_inversion
#' @export
ptweedie.inversion <- function(q, power, mu, phi, verbose, details){ 
  lifecycle::deprecate_warn(when = "3.0.5", 
                            what = "ptweedie.inversion()", 
                            with = "ptweedie_inversion()")
  ptweedie_inversion(q = q, 
                     power = power,
                     mu = mu, 
                     phi = phi, 
                     verbose = FALSE, 
                     details = FALSE)
}

