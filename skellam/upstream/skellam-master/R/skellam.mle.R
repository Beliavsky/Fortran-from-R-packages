#' Maximum Likelihood Estimation for the Skellam Distribution
#'
#' Estimates the parameters of a Skellam distribution using maximum likelihood.
#'
#' @param x A vector of integers (positive or negative).
#'
#' @details
#' Instead of having to maximize the log-likelihood with respect to both parameters
#' (\eqn{\lambda_1} and \eqn{\lambda_2}), the function maximizes with respect to
#' \eqn{\lambda_2} while setting \eqn{\lambda_1 = \lambda_2 + \bar{x}}. This
#' approach improves computational efficiency. The optimization is performed using
#' \code{\link[stats]{nlm}} as it proved faster than \code{\link[stats]{optimise}}.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{iters}{Number of iterations required by \code{nlm}.}
#'   \item{loglik}{Maximized log-likelihood value.}
#'   \item{param}{Estimated parameters (\eqn{\hat{\lambda}_1}, \eqn{\hat{\lambda}_2}).}
#' }
#'
#' @references
#' \itemize{
#'   \item Butler, R. (2007) \emph{Saddlepoint Approximations with Applications},
#'   Cambridge University Press.
#'
#'   \item Johnson, N. L. (1959) On an extension of the connection between Poisson
#'   and \eqn{\chi^2} distributions. \emph{Biometrika} \bold{46}, 352-362.
#'
#'   \item Johnson, N. L.; Kotz, S.; Kemp, A. W. (1993) \emph{Univariate Discrete
#'   Distributions}, 2nd ed., John Wiley and Sons.
#'
#'   \item Skellam, J. G. (1946) The frequency distribution of the difference between
#'   two Poisson variates belonging to different populations. \emph{Journal of the
#'   Royal Statistical Society, Series A} \bold{109}(3), 296.
#'
#'   \item Strackee, J.; van der Gon, J. J. D. (1962) The frequency distribution of
#'   the difference between two Poisson variates. \emph{Statistica Neerlandica}
#'   \bold{16}(1), 17-23.
#'
#'   \item Abdulhamid, A. A.; Maha, A. O. (2010) On The Poisson Difference
#'   Distribution Inference and Applications. \emph{Bulletin of the Malaysian
#'   Mathematical Sciences Society} \bold{33}(1), 17-45.
#'
#'   \item Wikipedia: Skellam distribution \url{https://en.wikipedia.org/wiki/Skellam_distribution}
#' }
#'
#' @examples
#' # Basic example
#' x1 <- rpois(1000, 10)
#' x2 <- rpois(1000, 6)
#' x <- x1 - x2
#' skellam.mle(x)
#'
#' # Larger sample size
#' x1 <- rpois(10000, 10)
#' x2 <- rpois(10000, 6)
#' x <- x1 - x2
#' skellam.mle(x)
#'
#' @export
#' @keywords distribution models
#' @author Michail Tsagris
#' @export
skellam.mle <- function(x) {

  n <- length(x)
  sx <- sum(x)
  sx2 <- sx / 2
  mx <- sx / n
  theta <- stats::var(x)/2  - mx/2

  skel <- function(theta) {
    a <- 2 * sqrt(theta^2 + theta * mx)
    n * 2 * theta + sx - sx2 * log(1 + mx / theta) -
    sum( log( besselI(a, x, expon.scaled = TRUE) ) ) - n * a
  }

  options(warn = -1)
  mod <- stats::nlm(skel, theta)
  param <- c( mod$estimate + mx, mod$estimate )
  names(param) <- c("mu 1", "mu 2")
  list(iters = mod$iterations, loglik = -mod$minimum, param = param)

}






