#'Rescaling Derivative Estimates
#'
#'The estimation functions of the \code{smoots} package estimate the
#'nonparametric trend function or its derivatives on the rescaled
#'time interval \eqn{[0, 1]}. With this function the derivative estimates can
#'be rescaled in accordance with a given vector with time points.
#'
#'@param y a numeric vector or matrix with the derivative estimates obtained
#'for time points on the interval \eqn{[0, 1]}; pass the list element \code{ye}
#'of the output of the functions \code{\link{dsmooth}} or \code{\link{gsmooth}}
#'(if the argument \code{v} \eqn{> 0}) to this argument.
#'@param x a numeric vector of length \code{length(y)} with the actual
#'(equidistant) time points ordered from past to present; the default is
#'\code{seq_along(y)}.
#'@param v the order of derivative that is implemented for \code{y}; the
#'default is \code{1}.
#'
#'@export
#'
#'@details
#'The derivative estimation process is based on the additive time series model
#'\deqn{y_t = m(x_t) + \epsilon_t,}
#'where \eqn{y_t} is the observed time series with equidistant design,
#'\eqn{x_t} is the rescaled time on \eqn{[0, 1]}, \eqn{m(x_t)} is a smooth and
#'deterministic trend function and \eqn{\epsilon_t} are stationary errors
#'with E(eps_[t]) = 0 (see also Beran and Feng, 2002). Since the estimates of
#'the main smoothing functions in \code{smoots} are obtained with regard to the
#'rescaled time points \eqn{x_t}, the derivative estimates returned by these
#'functions are valid for \eqn{x_t} only. Thus, by passing the returned
#'estimates to the argument \code{y}, a vector with the actual time points to
#'the argument \code{x} and the order of derivative of \code{y} to the argument
#'\code{v}, a rescaled estimate series is calculated and returned. The function
#'can also be combined with the numeric output of \code{\link{confBounds}}.
#'
#'Note that the trend estimates, even though they are also obtained for the
#'rescaled time points \eqn{x_t}, are still valid for the actual time points.
#'
#'@return
#'A numeric vector with the rescaled derivative estimates is returned.
#'
#'@author
#'\itemize{
#'\item Dominik Schulz (Research Assistant) (Department of Economics, Paderborn
#'University), \cr
#'Package Creator and Maintainer
#'}
#'
#'@examples
#'data <- smoots::gdpUS
#'Xt <- log(data$GDP)
#'time <- seq(from = 1947.25, to = 2019.5, by = 0.25)
#'d_est <- smoots::dsmooth(Xt)
#'ye_rescale <- smoots::rescale(d_est$ye, x = time, v = 1)
#'plot(time, ye_rescale, type = "l", main = "", ylab = "", xlab = "Year")
#'

rescale <- function(y, x = seq_along(y), v = 1) {

  x.diff <- x[2] - x[1]
  n <- length(x)
  scl <- x[n] - x[1] + x.diff
  y / (scl^v)
}
