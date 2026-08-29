#' @title Tweedie Distribution: Convert Between Parameter Formats
#' @name tweedie_convert
#' @description
#' Converts from the fitted \acronym{glm} parameters \eqn{p}, \eqn{\mu}{mu} and \eqn{\phi}{phi}
#' and the corresponding underlying Poisson and gamma parameters (when \eqn{1 < p < 2}).
#'
#' @param xi a synonym for \code{power}.
#' @param mu the mean parameter \eqn{\mu}{mu}.
#' @param phi the dispersion parameter \eqn{\phi}{phi}.
#' @param power the power parameter \eqn{p}{power}; a synonym for \eqn{\xi}{xi}.
#' @return a list of the parameters of the parameters of the corresponding underlying 
#' Poisson and gamma densities: 
#' \code{poisson.lambda} (\eqn{\lambda}{lambda} from the underlying Poisson distribution),
#' \code{gamma.shape}, \code{gamma.scale} (the shape and scale parameters
#'  from the underlying gamma distribution),
#'  \code{p0} (the probability that \eqn{Y = 0}),
#'  \code{gamma.mean} and \code{gamma.phi} (the gamma mean and dispersion parameter values)
#'  
#' @importFrom stats glm fitted
#' @importFrom statmod tweedie
#' 
#' @examples
#' ### Fit a Tweedie density
#' pretend <- data.frame( y = rgamma(20, shape = 1, rate = 1) )
#' fit <- glm(y ~ 1, data = pretend, 
#'            family = statmod::tweedie(link.power = 0, var.power = 1.4))
#'
#' # Convert parameters
#' tweedie_convert(mu = fitted(fit, type="response"), phi = 1, power = 1.4)
#'
#' @export
tweedie_convert <- function(xi = NULL, mu, phi, power = NULL){
  ### ADDED 14 July 2017
  
  if ( is.null(power) & is.null(xi) ) stop("Either xi or power must be given\n")
  xi.notation <- TRUE
  if ( is.null(power) ) {   # Then  xi  is given
    if ( !is.numeric(xi)) stop("xi  must be numeric.\n")
    power <- xi
  } else {
    xi.notation <- FALSE
  }
  if ( is.null(xi) ) {   # Then   power  is given
    if ( !is.numeric(power)) stop("power  must be numeric.\n")
    xi.notation <- FALSE
    xi <- power
  }
  if ( xi != power ) {
    cat("Different values for xi and power given; the value of xi used.\n")
    power <- xi
  }
  index.par       <- ifelse( xi.notation, "xi", "p")
  index.par.long  <- ifelse( xi.notation, "xi", "power")
  
  
  # Error checks
  if ( power < 1)      stop( paste(index.par.long, "must be greater than 1.\n") )
  if ( power >= 2)     stop( paste(index.par.long, "must be less than 2.\n") )
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(mu <= 0) )  stop("mu must be positive.\n")
  
  if( length(mu) != length(phi) ){
    if ( length(mu)  == 1 ) mu  <- array(dim = length(phi), mu  ) 
    if ( length(phi) == 1 ) phi <- array(dim = length(mu),  phi ) 
  }
  # Now  mu  and  phi  will be the same length if one of them happened to be a scalar, so this works:
  if( length(mu) != length(phi) ) stop("phi and mu must be scalars, or the same length.\n")
  
  lambda <- ( mu ^ (2 - xi) ) / ( phi * (2 - xi) )  # Poisson distribution mean
  alpha  <- (2 - xi)  / (xi - 1)             # gamma distribution alpha (shape)
  gam    <- phi * (xi - 1) * mu ^ (xi - 1)   # gamma distribution beta  (scale)
  p0     <- exp( -lambda )
  phi.g  <- (2 - xi) * (xi - 1) * phi ^ 2 * mu ^ ( 2 * (xi-1) )
  mu     <- gam / phi
  
  list( poisson.lambda = lambda, 
        gamma.shape = alpha, 
        gamma.scale = gam, 
        p0 = p0,
        gamma.mean = mu,
        gamma.phi = phi.g)
  
}





#' @rdname tweedie_convert
#' @export
tweedie.convert <- function(xi = NULL, mu, phi, power = NULL){
  lifecycle::deprecate_warn(when = "3.0.5", 
                            what = "tweedie.convert()", 
                            with = "tweedie_convert()")
  if (is.null(power)) power <- xi
  tweedie_convert(xi = NULL, 
              mu = mu,
              phi = phi,
              power=power)
}
