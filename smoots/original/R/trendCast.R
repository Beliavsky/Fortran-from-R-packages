#' Forecasting Function for Nonparametric Trend Functions
#'
#'@param object an object returned by either \code{\link{msmooth}},
#'\code{\link{tsmooth}}, \code{\link{gsmooth}} (with \code{v = 0}) or
#'\code{\link{knsmooth}}.
#'@param h the forecasting horizon; the values \eqn{m(n + 1)} to \eqn{m(n + h)}
#'will be predicted; is set to \code{h = 1} by default; decimal numbers will be
#'rounded off to integers.
#'@param np.fcast the forecasting method; \code{np.fcast = "lin"} uses a linear
#'extrapolation, whereas \code{np.fcast = "const"} uses the last fitted value
#'of \eqn{m(x_t)} as a forecast; is set to \code{"lin"} by default.
#'@param plot a logical value; if set to \code{TRUE}, a simple plot of the
#'original time series, the local polynomial trend estimates as well as the
#'predicted values is generated.
#'@param ... additional arguments for the standard plot function, e.g.,
#'\code{xlim}, \code{type}, ... ; arguments with respect to plotted graphs,
#'e.g., the argument \code{col}, only affect the original series \code{X};
#'please note that in accordance with the argument \code{x} (lower case) of the
#'standard plot function, an additional numeric vector with time points can be
#'implemented via the argument \code{x} (lower case). \code{x} should be
#'valid for the sample observations only, i.e.
#'\code{length(x) == length(obj$orig)} should be \code{TRUE}, as future time
#'points will be calculated automatically.
#'
#'@details
#'This function is part of the \code{smoots} package and was implemented under
#'version 1.1.0. The underlying theory is based on the additive nonparametric
#'regression function
#'\deqn{y_t = m(x_t) + \epsilon_t,}
#'where \eqn{y_t} is the observed time series with equidistant design,
#'\eqn{x_t} is the rescaled time on the interval \eqn{[0, 1]}, \eqn{m(x_t)}
#'is a smooth and deterministic trend function and \eqn{\epsilon_t} are
#'stationary errors with \eqn{E(\epsilon_t) = 0}.
#'
#'The purpose of this function is the forecasting of future values based on
#'a nonparametric regression model. Following the proposition in Fritz
#'et al. (2020), point predictions can be conducted
#'separately for the nonparametric trend function \eqn{m(x_t)} and the
#'stationary part \eqn{\epsilon_t}. The sum of both forecasts is then the
#'forecast of \eqn{y_t}. With this function, only the forecast with respect to
#'\eqn{m(x_t)} is computable.
#'Now assume that the variance of the error in the local polynomial
#'forecasts is negligible when calculating the forecasting intervals. We
#'define the forecast for time point \eqn{n + k}, \eqn{k = 1, 2, ..., h}, by
#'\deqn{\hat{m}(x_{n + k}) = \hat{m}(x_n) + D k \delta_m,}{hat[m](x_[n + k]) =
#'hat[m](x_n) + Dk * \delta_m,}
#'where \eqn{\delta_m} is equal to \eqn{\hat{m}(x_n) -
#'\hat{m}(x_{n - 1})}{hat[m](x_n) - hat[m](x_[n - 1])} and \eqn{D} is a
#'dummy  variable. If \eqn{D = 1}, a linear extrapolation is applied. For
#'\eqn{D = 0}, \eqn{\hat{m}(x_n)}{hat[m](x_n)} is the predicted value.
#'
#'To make use of this function, an object of class \code{smoots} can be given
#'as input. However, since the discussed approach is only valid for the
#'estimated trend function, only objects created by \code{\link{msmooth}},
#'\code{\link{tsmooth}}, \code{\link{knsmooth}} and \code{link{gsmooth}}, if
#'the trend was estimated, will be appropriate input objects.
#'
#'With the input argument \code{h}, a positive integer can be given to the
#'function that represents the forecasting horizon, i.e. how many future values
#'are to be estimated. Via the argument \code{np.fcast} the value of the dummy
#'variable D can be specified and thus the forecasting method. For
#'\code{np.fcast = "lin"}, \eqn{D = 1} is applied, whereas for
#'\code{np.fcast = "const"}, \eqn{D} is set to \eqn{0}.
#'
#'By means of the argument \code{plot} that can be either set to the logical
#'values \code{TRUE} or \code{FALSE}, a simple plot of the original series
#'alongside the local polynomial estimates as well as the forecasted values can
#'be either generated or suppressed.
#'
#'The function always returns a vector of forecasted values ordered from
#'\eqn{n + 1} to \eqn{n + h}. Depending on the setting of the argument
#'\code{plot}, a generic plot of the results may be generated. Furthermore,
#'additional arguments of the standard plot function can be passed to this
#'function as well to adjust the generated plot.
#'
#'@export
#'
#'@importFrom graphics lines
#'@importFrom utils tail
#'
#'@return
#'A numeric vector is always returned with the forecasted values. Depending
#'on the setting for the argument \emph{plot}, a generic plot might be created.
#'
#'@references
#' Feng, Y., Gries, T. and Fritz, M. (2020). Data-driven
#' local polynomial for the trend and its derivatives in economic time
#' series. Journal of Nonparametric Statistics, 32:2, 510-533.
#'
#' Feng, Y., Gries, T., Letmathe, S. and Schulz, D. (2019). The smoots package
#' in R for semiparametric modeling of trend stationary time series. Discussion
#' Paper. Paderborn University. Unpublished.
#'
#' Feng, Y., Gries, T., Fritz, M., Letmathe, S. and Schulz, D. (2020).
#' Diagnosing the trend and bootstrapping the forecasting intervals using a
#' semiparametric ARMA. Discussion Paper. Paderborn University. Unpublished.
#'
#' Fritz, M., Forstinger, S., Feng, Y., and Gries, T. (2020). Forecasting
#' economic growth processes for developing economies. Unpublished.
#'
#'@author
#'\itemize{
#'\item Yuanhua Feng (Department of Economics, Paderborn University), \cr
#'Author of the Algorithms \cr
#'Website: \url{https://wiwi.uni-paderborn.de/en/dep4/feng/}
#'\item Dominik Schulz (Research Assistant) (Department of Economics, Paderborn
#'University), \cr
#'Package Creator and Maintainer
#'}
#'
#'@examples
#'log_gdp <- log(smoots::gdpUS$GDP)
#'est <- msmooth(log_gdp)
#'forecasts <- trendCast(est, h = 5, plot = TRUE)
#'forecasts
#'

trendCast <- function(object, h = 1, np.fcast = c("lin", "const"),
  plot = FALSE, ...) {

  if (!inherits(object, "smoots")) {
    stop("Input object not of class 'smoots'.")
  }
  if (attr(object, "function") %in% c("dsmooth", "confBounds")) {
    stop("Only output objects of trend estimation functions are allowed as ",
      "input.")
  }
  if (attr(object, "function") == "gsmooth" && object$v != 0) {
    stop("Only output objects of trend estimation functions are allowed as ",
      "input.")
  }
  if (length(h) != 1 || is.na(h) || !is.numeric(h) || h < 1) {
    stop("Argument 'h' must be a positive integer or at least a value >= 1.")
  }
  h <- floor(h)

  if (!(length(np.fcast) %in% c(1, 2)) || !all(!is.na(np.fcast)) ||
      !is.character(np.fcast)) {
    stop("Argument 'np.fcast' is defined incorrectly.")
  }
  if (all(np.fcast == c("lin", "const"))) np.fcast <- "lin"
  if (length(np.fcast) != 1 || !(np.fcast %in% c("lin", "const"))) {
    stop("Argument 'np.fcast' is defined incorrectly.")
  }

  if (length(plot) != 1 || is.na(plot) || !(plot %in% c(TRUE, FALSE))) {
    stop("Argument 'plot' must be a logical value (TRUE or FALSE).")
  }

  ye <- object$ye
  y <- object$orig
  n <- object$n

  if (np.fcast == "lin") {
    b <- ye[n] - ye[n - 1]
    out <- ye[n] + (1:h) * b
  } else {
    out <- rep(ye[n], times = h)
  }

  if (plot == TRUE) {
    dots <- list(...)
    if (is.null(dots[["xlab"]])) dots[["xlab"]] <- "Time"
    if (is.null(dots[["ylab"]])) dots[["ylab"]] <- "Series"
    if (is.null(dots[["main"]])) dots[["main"]] <- "Forecasting results"
    if (is.null(dots[["type"]])) dots[["type"]] <- "l"
    if (is.null(dots[["x"]])) dots[["x"]] <- 1:n
    if (is.null(dots[["y"]])) dots[["y"]] <- y
    xStep <- dots$x[2] - dots$x[1]
    x.fc <- seq(dots$x[[n]] + xStep, dots$x[[n]] + h * xStep, xStep)
    if (is.null(dots[["xlim"]])) {
      xlim1 <- tail(dots$x, 6 * h)[1]
      dots[["xlim"]] <- c(xlim1, dots$x[[n]] + h * xStep)
    }
    if (is.null(dots[["ylim"]])) {
      x.all <- c(dots$x, x.fc)
      bound.y <- c(y, out)
      bound.ye <- c(ye, out)
      x.low <- sum(x.all < dots$xlim[1]) + 1
      x.up <- sum(x.all <= dots$xlim[2])
      dots[["ylim"]] <- c(min(bound.y[x.low:x.up], bound.ye[x.low:x.up]),
        max(bound.y[x.low:x.up], bound.ye[x.low:x.up]))
    }

    do.call(what = graphics::plot, args = dots)
    lines(dots$x, ye, col = "red", lty = 1)
    lines(c(dots$x[[n]], x.fc), c(ye[n], out), col = "deepskyblue4", lty = 1)
  }
  out
}
