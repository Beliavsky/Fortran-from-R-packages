#' @title Fourier Inversion Evaluation for the Tweedie Probability Function
#' @name dtweedie_inversion
#' @description
#' Evaluates the probability density function (\acronym{pdf}) for Tweedie distributions using Fourier inversion, 
#' for given values of the dependent variable \code{y}, the mean \code{mu}, dispersion \code{phi}, and power parameter \code{power}.
#' \emph{Not usually called by general users}, but can be used in the case of evaluation problems.
#'
#' @usage dtweedie_inversion(y, mu, phi, power, method = 3, verbose = FALSE, 
#'                           details = FALSE, IGexact = TRUE)
#'
#' @param y vector of quantiles.
#' @param mu the mean parameter \eqn{\mu}{mu}.
#' @param phi the dispersion parameter \eqn{\phi}{phi}.
#' @param power scalar; the power parameter \eqn{p}{power}.
#' @param method the method to use; one of \code{1}, \code{2}, or \code{3} (the default).
#' @param verbose logical; if \code{TRUE}, display some internal computation details. The default is \code{FALSE}.
#' @param details logical; if \code{TRUE}, return a list with basic details of the integration. The default is \code{FALSE}.
#' @param IGexact logical; if \code{TRUE} (the default), evaluate the inverse Gaussian distribution using the 'exact' values, otherwise uses inversion.
#'
#' @return A numeric vector of densities if \code{details=FALSE}; if \code{details = TRUE}, a list containing \code{denisty} (a vector of the values of the density), \code{regions} (a vector of the number of integration regions used),\code{method} (a vector giving the evaluation method used; see the Note below on the three methods), and \code{exitstatus} (a vector, where a \code{1} for any value means a computational problem or target relative accuracy not reached, for the corresponding observation).
#' 
#' @note
#' The 'exact' values for the inverse Gaussian distribution are not really exact, but evaluated using inverse normal distributions,
#' for which very good numerical approximation are available in R.

#' For special cases of \eqn{p} (i.e., \eqn{p = 0, 1, 2, 3}), where no inversion is needed, \code{regions} and \code{method} are set to \code{NA} for all values of \code{y}.
#' For special cases of \code{y} for other values of \eqn{p} (i.e., \eqn{P(Y = 0)}), \code{regions} and \code{method} are set to \code{NA}.
#'
#' @note
#' The three methods are described in Dunn & Smyth (2008).
#' 
#' @references
#' Dunn, P. K. and Smyth, G. K. (2008).
#' Evaluation of Tweedie exponential dispersion model densities by Fourier inversion.
#' \emph{Statistics and Computing}, 
#' \bold{18}, 73--86.
#' \doi{10.1007/s11222-007-9039-6}

#' @examples
#' # Plot a Tweedie density
#' y <- seq(0.02, 4, length = 50)
#' fy <- dtweedie_inversion(y, mu = 1, phi = 1, power = 1.1)
#' plot(y, fy, type = "l", lwd = 2, ylab = "Density")
#' 
#' @references
#' Dunn, P. K. and Smyth, G. K. (2008).
#' Evaluation of Tweedie exponential dispersion model densities by Fourier inversion.
#' \emph{Statistics and Computing}, 
#' \bold{18}, 73--86.
#' \doi{10.1007/s11222-007-9039-6}
#' 
#' @keywords distribution
#' 
#' @export
dtweedie_inversion <- function(y, mu, phi, power, method = 3, verbose = FALSE,  
                               details = FALSE, IGexact = TRUE){ 
  ### NOTE: No notation checks
  
  # Check
  if (length(y) == 0L) {
    return(numeric(0))
  }
  
  # CHECK THE INPUTS ARE OK AND OF CORRECT LENGTHS
  if (verbose) cat("- Checking, resizing inputs\n")
  out <- check_inputs(y, mu, phi, power)
  mu <- out$mu
  phi <- out$phi
  
  # cdf    is the whole vector; the same length as  y.
  # All is resolved in the end.
  density <- numeric(length = length(y) )
  regions <- rep(NA, length(y)) 
  
  # IDENTIFY SPECIAL CASES
  special_y_cases <- rep(FALSE, length(y))
  if (verbose) cat("- Checking for special cases\n")
  out <- special_cases(y, mu, phi, power,
                       IGexact = IGexact,
                       type = "PDF")
  
  special_p_cases <- out$special_p_cases
  special_y_cases <- out$special_y_cases
  
  if (verbose & special_p_cases) cat("  - Special case for p used\n")
  if ( any(special_y_cases) ) {
    special_y_cases <- out$special_y_cases  
    if (verbose) cat("  - Special cases for first input found\n")
    density <- out$f # This is the final vector of results to return, filled with the special-case info.
    # NOTE: regions filled with zeros by default, so regions = 0 in these cases
  }
  
  if ( special_p_cases ) {
    density <- out$f
    optimal_Method <- array(NA, dim = length(y)) 
  } else {
		# NOT special p case; ONLY special y cases 
		
		# Now use FORTRAN on the remaining values:
		regions <- integer(length = length(y)) # Filled with zeros by default
		N <- length(y)
		N_nonSpecial <- N - sum(out$special_y_cases) 
		regions[special_y_cases] <- NA
    optimal_Method <- array(NA, dim = N) 

		
		### BEGIN SET UP
		pSmall  <- ifelse( (power > 1) & (power < 2), 
											 TRUE, FALSE )
	
		# Initialise
		exitstatus_scalar <- as.integer(0)
		relerr_scalar     <- as.double(0.0)
		its_scalar        <- as.integer(0)
		### END SET UP
		
		
		# Establish which method to use
		# For special cases, method is 0.
		if ( is.null(method)){
			method <- array(0, dim = N)
		} else {
			method <- array( method, 
											 dim = N)
		}
		method[special_y_cases] <- 0
		# There are three approaches ('method'), each a product of a simple bit
		# and a complicated bit computed in FORTRAN
		#
		# The methods are documented in Dunn and Smyth (2008):
		# - Method 1: Evaluate a(): compute a(y, phi) = f(y; 1, phi)
		# - Method 2: Rescale the mean to 1
		# - Method 3: Rescale y to 1 and evaluate b().
		#
		# If no method is explicitly requested, find the notional "optimal" method for each i.
		
		### BEGIN ESTABLISH METHOD
		theta <- ( mu ^ (1 - power) - 1 ) / ( 1 - power )
		if ( ( abs(power - 2 ) ) < 1.0e-07 ){
			kappa <- log(mu) + (2 - power) * ( log(mu) ^ 2 ) / 2
		} else {
			kappa <- ( mu ^ (2 - power) - 1 ) / (2 - power)
		}
	
		# Method 1
		m1 <- exp( (y * theta - kappa ) / phi )
		dev <- tweedie_dev(y = y, 
											 mu = mu,
											 power = power )
		
		# Method 2
		m2 <- 1 / mu
	
		# Method 3
		m3 <- exp( -dev/(2 * phi) ) / y
			
		# Select method: this is an n x 3 array of the vaklues of [m1, m2, m3], 
		# and from this we chose the method (i.e., column) containing the minimum
		method_List <- array(c(m1, m2, m3), 
												 dim = c(length(y), 3))
		method_List[special_y_cases, ] <- 0 # Method 0 for special-cases rows
	
		optimal_Method <- apply(method_List, 
														MARGIN = 1, 
														FUN = which.min)
		optimal_Method[special_y_cases] <- 0  
	
		### BEGIN: Set parameters for FORTRAN call, depending on method
		# mu = 1 for all methods:
		mu_F <- rep(1, length(y) )
		
		# Set up empty vector to fill for other methods:
		phi_F <- phi
		y_F <- y
	
		# Method 1 just uses the given  y  and  phi
		if (any(optimal_Method == 2)){
			use_M2 <- optimal_Method==2
			phi_F[ use_M2 ] <- phi[use_M2] / mu[use_M2] ^ (2 - power)
			y_F[ use_M2 ] <- y[use_M2]/mu[use_M2]
		}

		if (any(optimal_Method == 3)){
			use_M3 <- optimal_Method==3
			phi_F[ use_M3 ] <- phi[use_M3] / y[use_M3] ^ (2 - power)
			y_F[ use_M3 ] <- 1
		}
		### END: Set parameters for FORTRAN call, depending on method
	
		tmp <- .C("twcomputation",
							N          = as.integer(N_nonSpecial),
							power      = as.double(power),
							phi        = as.double(phi_F[!special_y_cases]),
							y          = as.double(y_F[!special_y_cases]),
							mu         = as.double(mu_F[!special_y_cases]),
							verbose    = as.integer( verbose ),
							pdf        = as.integer(1),          # 1: TRUE, as this is the PDF
							# THE OUTPUTS:
							funvalue   = numeric(N_nonSpecial),  # funvalue
							exitstatus = integer(N_nonSpecial),  # exitstatus
							relerr     = numeric(N_nonSpecial),  # relerr
							its        = integer(N_nonSpecial),  # its
							PACKAGE    = "tweedie")
		
		density[!special_y_cases] <- tmp$funvalue
		regions[!special_y_cases] <- tmp$its

		# Reconstruct
		if (any(optimal_Method == 1)){
			use_M1 <- optimal_Method==1
			density[use_M1] <- density[use_M1] * m1[use_M1]
		}    
		if (any(optimal_Method == 2)){
			density[use_M2] <- density[use_M2] * m2[use_M2]
		}  
		if (any(optimal_Method == 3)){
			density[use_M3] <- density[use_M3] * m3[use_M3]
		}
		
		# Now for special cases: set method to NA
		optimal_Method[special_y_cases] <- NA
  }

  # Return
  if (details) {
    return( list( density = density,
                  regions = regions,
                  method = optimal_Method,
                  exitstatus = tmp$exitstatus))
  } else {
    return(density)
  }
}


#' @rdname dtweedie_inversion
#' @export
dtweedie.inversion <- function(y, power, mu, phi, method = 3, verbose, details){ 
  lifecycle::deprecate_warn(when = "3.0.5", 
                            what = "dtweedie.inversion()", 
                            with = "dtweedie_inversion()")
  dtweedie_inversion(y = y, 
                     mu = mu, 
                     phi = phi, 
                     power = power, 
                     method = method, 
                     verbose = FALSE, 
                     details = FALSE)
}

