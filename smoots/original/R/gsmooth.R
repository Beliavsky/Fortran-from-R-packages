#' Estimation of Trends and their Derivatives via Local Polynomial Regression
#'
#' This function is an R function for estimating the trend function
#' and its derivatives in an equidistant time series with local polynomial
#' regression and a fixed bandwidth given beforehand.
#'
#' @param y a numeric vector that contains the time series data ordered from
#' past to present.
#' @param v an integer \code{0}, \code{1}, ... that represents the order of
#' derivative that will be estimated; is set to \code{v = 0} by default.
#'
#' \tabular{cl}{
#' \strong{Number (\code{v})} \tab \strong{Degree of derivative}\cr
#' \code{0} \tab The function \emph{f(x)} itself\cr
#' \code{1} \tab The first derivative \emph{f'(x)}\cr
#' \code{2} \tab The second derivative \emph{f''(x)}\cr
#' \code{...} \tab ...
#' }
#' @param p an integer \eqn{>= (} \code{v} \eqn{+ 1)} that represents the order
#' of polynomial; \code{p - v} must be an odd number; is set to \code{v + 1}
#' by default.
#'
#' Exemplary for \code{v = 0}:
#'
#' \tabular{clcll}{
#' \strong{Number (\code{p})} \tab \strong{Polynomial} \tab
#' \strong{\code{p - v}} \tab \strong{\code{p - v} odd?} \tab
#' \strong{\code{p} usable?}\cr
#' \code{1} \tab Linear \tab 1 \tab Yes
#' \tab Yes\cr
#' \code{2} \tab Quadratic \tab 2 \tab No
#' \tab No\cr
#' \code{3} \tab Cubic \tab 3 \tab Yes
#' \tab Yes\cr
#' \code{...} \tab ... \tab ... \tab ...
#' \tab ...
#' }
#' @param mu an integer \code{0}, \code{1}, \code{2}, ... that represents the
#' smoothness parameter of the kernel weighting function that will be used; is
#' set to \code{1} by default.
#'
#' \tabular{cl}{
#'\strong{Number (\code{mu})} \tab \strong{Kernel}\cr
#'\code{0} \tab Uniform Kernel\cr
#'\code{1} \tab Epanechnikov Kernel\cr
#'\code{2} \tab Bisquare Kernel\cr
#'\code{3} \tab Triweight Kernel\cr
#'\code{...} \tab ...
#' }
#' @param b a real number \eqn{0 <} \code{b} \eqn{< 0.5}; represents the
#' relative bandwidth that will be used for the smoothing process; is set to
#' \code{0.15} by default.
#' @param bb can be set to \code{0} or \code{1}; the parameter controlling the
#' bandwidth used at the boundary; is set to \code{1} by default.
#'
#' \tabular{cl}{
#' \strong{Number (\code{bb})} \tab \strong{Estimation procedure at boundary
#' points}\cr
#' \code{0} \tab Fixed bandwidth on one side with possible large
#' bandwidth on the other side at the boundary\cr
#' \code{1} \tab The k-nearest neighbor method will be used
#' }
#'
#'@export
#'
#'@details
#'The trend or its derivatives are estimated based on the additive
#'nonparametric regression model for an equidistant time series
#'\deqn{y_t = m(x_t) + \epsilon_t,}
#'where \eqn{y_t} is the observed time series, \eqn{x_t} is the rescaled time
#'on the interval \eqn{[0, 1]}, \eqn{m(x_t)} is a smooth and deterministic
#'trend function and \eqn{\epsilon_t} are stationary errors with
#'\eqn{E(\epsilon_t) = 0} (see also Beran and Feng, 2002).
#'
#'This function is part of the package \code{smoots} and is used in
#'the field of analyzing equidistant time series data. It applies the local
#'polynomial regression method to the input data with an arbitrarily
#'selectable bandwidth. By these means, the trend as well as its derivatives
#'can be estimated nonparametrically, even though the result will strongly
#'depend on the bandwidth given beforehand as an input.
#'
#'NOTE:
#'
#'The estimates are obtained with regard to the rescaled time points on the
#'interval \eqn{[0, 1]}. Thus, if \eqn{\nu > 0}, the estimates might not
#'reflect the values for the actual time points. To rescale the estimates, we
#'refer the user to the \code{\link{rescale}} function of the \code{smoots}
#'package.
#'
#'With package version 1.1.0, this function implements C++ code by means
#'of the \code{\link[Rcpp:Rcpp-package]{Rcpp}} and
#'\code{\link[RcppArmadillo:RcppArmadillo-package]{RcppArmadillo}} packages for
#'better performance.
#'
#'@md
#'
#'@return The output object is a list with different components:
#'
#'\describe{
#'\item{b}{the chosen (relative) bandwidth; input argument.}
#'\item{bb}{the chosen bandwidth option at the boundaries; input argument.}
#'\item{mu}{the chosen smoothness parameter for the second order kernel; input
#'argument.}
#'\item{n}{the number of observations.}
#'\item{orig}{the original input series; input argument.}
#'\item{p}{the chosen order of polynomial; input argument.}
#'\item{res}{a vector with the estimated residual series; is set to \code{NULL}
#'for \code{v > 0}.}
#'\item{v}{the order of derivative; input argument.}
#'\item{ws}{the weighting system matrix used within the local polynomial
#'regression; this matrix is a condensed version of a complete weighting system
#'matrix; in each row of \code{ws}, the weights for conducting the smoothing
#'procedure at a specific observation time point can be found; the first
#'\eqn{[nb + 0.5]} rows, where \eqn{n} corresponds to the number of
#'observations, \eqn{b} is the bandwidth considered for smoothing and
#'\eqn{[.]} denotes the integer part, contain the weights at the
#'\eqn{[nb + 0.5]} left-hand boundary points; the weights in row
#'\eqn{[nb + 0.5] + 1} are representative for the estimation at all
#'interior points and the remaining rows contain the weights for the right-hand
#'boundary points; each row has exactly \eqn{2[nb + 0.5] + 1} elements,
#'more specifically the weights for observations of the nearest
#'\eqn{2[nb + 0.5] + 1} time points; moreover, the weights are normalized,
#'i.e. the weights are obtained under consideration of the time points
#'\eqn{x_t = t/n}, where \eqn{t = 1, 2, ..., n}.}
#'\item{ye}{a vector with the estimates of the selected nonparametric order of
#'derivative on the rescaled time interval \eqn{[0, 1]}.}
#'}
#'
#'@references
#' Beran, J. and Feng, Y. (2002). Local polynomial fitting with long-memory,
#' short-memory and antipersistent errors. Annals of the Institute of
#' Statistical Mathematics, 54(2), 291-311.
#'
#' Feng, Y., Gries, T. and Fritz, M. (2020). Data-driven
#' local polynomial for the trend and its derivatives in economic time
#' series. Journal of Nonparametric Statistics, 32:2, 510-533.
#'
#' Feng, Y., Gries, T., Letmathe, S. and Schulz, D. (2019). The smoots package
#' in R for semiparametric modeling of trend stationary time series. Discussion
#' Paper. Paderborn University. Unpublished.
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
#'# Logarithm of test data
#'test_data <- gdpUS
#'y <- log(test_data$GDP)
#'
#'# Applied gsmooth function for the trend with two different bandwidths
#'results1 <- gsmooth(y, v = 0, p = 1, mu = 1, b = 0.28, bb = 1)
#'results2 <- gsmooth(y, v = 0, p = 1, mu = 1, b = 0.11, bb = 1)
#'trend1 <- results1$ye
#'trend2 <- results2$ye
#'
#'# Plot of the results
#'t <- seq(from = 1947, to = 2019.25, by = 0.25)
#'plot(t, y, type = "l", xlab = "Year", ylab = "log(US-GDP)", bty = "n",
#'  lwd = 2,
#'  main = "Estimated trend for log-quarterly US-GDP, Q1 1947 - Q2 2019")
#'points(t, trend1, type = "l", col = "red", lwd = 1)
#'points(t, trend2, type = "l", col = "blue", lwd = 1)
#'legend("bottomright", legend = c("Trend (b = 0.28)", "Trend (b = 0.11)"),
#'  fill = c("red", "blue"), cex = 0.6)
#'title(sub = expression(italic("Figure 1")), col.sub = "gray47",
#'  cex.sub = 0.6, adj = 0)
#'
#'

# The basic function "gsmooth(y,v,p,mu,b,bb)"-------------------------------

gsmooth <- function(y, v = 0, p = v + 1, mu = 1, b = 0.15, bb = c(0, 1)) {

    if (length(y) <= 1 || !all(!is.na(y)) || !is.numeric(y)) {
      stop("The argument 'y' must be a numeric vector with length > 1 and ",
           "without NAs.")
    }

    if (length(v) != 1 || is.na(v) || !is.numeric(v) || v < 0) {
      stop("The argument 'v' must be a single non-zero and non-NA integer ",
           "value.")
    }
    v <- floor(v)

    if (length(p) != 1 || is.na(p) || !is.numeric(p) || p < 1) {
      stop("The argument 'p' must be a single strictly positive non-NA ",
        "integer value.")
    }
    p <- floor(p)

    if (v >= p || ((p - v) %% 2) == 0) {
      stop("p must be > v and p - v must be odd.")
    }

    if (length(mu) != 1 || is.na(mu) || !is.numeric(mu) || mu < 0) {
      stop("The argument 'mu' must be a single non-zero and non-NA integer ",
           "value.")
    }
    mu <- floor(mu)

    if (length(b) != 1 || is.na(b) || !is.numeric(b) || (b <= 0 || b >= 0.5)) {
      stop("The argument 'b' must be a single non-NA double value with ",
           "0 < b < 0.5.")
    }

    if (!length(bb) %in% c(1, 2) || !all(!is.na(bb)) || !is.numeric(b)) {
      stop("The argument 'bb' must be a single non-NA integer ",
           "(either 0 or 1).")
    }
    bb <- floor(bb)
    if (all(bb == c(0, 1))) bb <- 1
    if (length(bb) != 1 || !(bb %in% c(0, 1))) {
      stop("The argument 'bb' must be a single non-NA integer ",
           "(either 0 or 1).")
    }

    estim <- gsmoothCalc2Cpp(y, v, p, mu, b, bb)
    gr <- c(estim$ye)
    ws <- estim$ws
    result <- list(ye = gr, orig = y, ws = ws, v = v, p = p, mu = mu, b = b,
      bb = bb, n = length(y))
    if (v == 0) {
      result[["res"]] <- y - gr
    } else {
      result[["res"]] <- NULL
    }

    class(result) <- "smoots"  # set package attribute
    attr(result, "function") <- "gsmooth"  # set function attribute for
                                                 # own print function

    result
}
# End of the function
