#' @title Unit Deviance for a Tweedie Distribution
#' @name tweedie_dev
#' @description Computes the unit deviance for Tweedie distributions.
#'
#' @usage tweedie_dev(y, mu, power)
#' @param y vector of quantiles.
#' @param power the power parameter \eqn{p}{power}.
#' @param mu the mean parameter \eqn{\mu}{mu}.
#' 
#' @return A numeric vector containing the unit deviance.
#' 
#' @references
#' Jorgensen, B. (1997).
#' \emph{Theory of Dispersion Models}.
#' Chapman and Hall, London.
#' 
#' @examples
#' # Unit deviance is not symmetric in general:
#' round( tweedie_dev(0:6, mu = 3, power = 1.1), 3)
#' 
#' @export
tweedie_dev <- function(y, mu, power){
  # 
  # Peter K Dunn 
  # 29 Oct 2000 
  # 
  
  p <- power
  if(p == 1) {
    dev <- array( dim = length(y) )
    mu <- array( dim = length(y), mu )
    dev[y != 0] <- y[y != 0] * log(y[y != 0] / mu[y != 0]) - (y[y != 0] - mu[y != 0])
    dev[y == 0] <- mu[y == 0]
  } else{
    if(p == 2) {
      dev <- log(mu / y) + (y / mu) - 1
    } else{
      if (p == 0) {
        dev <- (y - mu) ^ 2
        dev <- dev / 2
      } else{
        dev <- (y ^ (2 - p)) / ((1 - p) * (2 - p)) - 
          (y * (mu ^ (1 - p)))/(1 - p) + 
          (mu ^ (2 - p)) / (2 - p)
      }
    }
  }
  dev * 2
}



#' @rdname tweedie_dev
#' @export
tweedie.dev <- function(y, mu, power){ 
  lifecycle::deprecate_warn(when = "3.0.5", 
                            what = "tweedie.dev()", 
                            with = "tweedie_dev()")
  tweedie_dev(y, mu, power)
}

