#' Mean and variance for hyper-Poisson distribution
#'
#' @author Freddy Hernandez, \email{fhernanb@unal.edu.co}
#'
#' @description
#' This function calculates the mean and variance for the
#' hyper-Poisson distribution with parameters \eqn{\mu} and \eqn{\sigma}.
#'
#' @param mu value of the mu parameter.
#' @param sigma value of the sigma parameter.
#'
#' @references
#' Sáez-Castillo, A. J., & Conde-Sánchez, A. (2013). A hyper-Poisson regression model for overdispersed and underdispersed count data. Computational Statistics & Data Analysis, 61, 148-157.
#'
#' @seealso \link{HYPERPO}.
#'
#' @details
#' The hyper-Poisson distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and density given by
#'
#' \eqn{f(x | \mu, \sigma) = \frac{\mu^x}{_1F_1(1;\mu;\sigma)}\frac{\Gamma(\sigma)}{\Gamma(x+\sigma)}}
#'
#' where the function \eqn{_1F_1(a;c;z)} is defined as
#'
#' \eqn{_1F_1(a;c;z) = \sum_{r=0}^{\infty}\frac{(a)_r}{(c)_r}\frac{z^r}{r!}}
#'
#' and \eqn{(a)_r = \frac{\gamma(a+r)}{\gamma(a)}} for \eqn{a>0} and \eqn{r} positive integer.
#'
#' This function calculates the mean and variance of this distribution.
#'
#' @return
#' the function returns a list with the mean and variance.
#'
#' @example examples/examples_mean_var_hp.R
#'
#' @name mean_var_hp
#'
#' @rdname mean_var_hp
#' @export
mean_var_hp <- function(mu, sigma) {
  media <- mu - (sigma - 1) * (F11(c=sigma, z=mu)-1) / F11(c=sigma, z=mu)
  varia <- mu + (mu-(sigma-1)) * media - media^2
  return(list(mean=media, variance=varia))
}
#'
#' @rdname mean_var_hp
#' @export
mean_var_hp2 <- function(mu, sigma) {
  media <- mu
  lambda <- obtaining_lambda(media=mu, gamma=sigma)
  varia <- lambda + (lambda-(sigma-1)) * media - media^2
  return(list(mean=media, variance=varia))
}

