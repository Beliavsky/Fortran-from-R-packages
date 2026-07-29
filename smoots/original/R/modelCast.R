#' Forecasting Function for Trend-Stationary Time Series
#'
#'Point forecasts and the respective forecasting intervals for
#'trend-stationary time series are calculated.
#'
#'@param obj an object of class \code{smoots}; must be the output of a trend
#'estimation process and not of a first or second derivative estimation
#'process.
#'@param p an integer value \eqn{>= 0} that defines the AR order \eqn{p} of the
#'underlying ARMA(\eqn{p,q}) model within the rest term (see the section
#'\emph{Details} for more information); is set to \code{NULL} by default; if no
#'value is passed to \code{p} but one is passed to \code{q}, \code{p} is set to
#'\code{0}; if both \code{p} and \code{q} are \code{NULL}, optimal orders
#'following the BIC for \eqn{0 \leq p,q \leq 5}{0 \le p,q \le 5} are chosen; is
#'set to \code{NULL} by default; decimal numbers will be rounded off to
#'integers.
#'@param q an integer value \eqn{\geq 0}{\ge 0} that defines the MA order
#'\eqn{q} of the underlying ARMA(\eqn{p,q}) model within \code{X}; is set to
#'\code{NULL} by default; if no value is passed to \code{q} but one is passed
#'to \code{p}, \code{q} is set to \code{0}; if both \code{p} and \code{q} are
#'\code{NULL}, optimal orders following the BIC for
#'\eqn{0 \leq p,q \leq 5}{0 \le p,q \le 5} are chosen; is set to \code{NULL} by
#'default; decimal numbers will be rounded off to integers.
#'@param h an integer that represents the forecasting horizon; if \eqn{n} is
#'the number of observations, point forecasts and forecasting intervals will be
#'obtained for the time points \eqn{n + 1} to \eqn{n + h}; is set to
#'\code{h = 1} by default; decimal numbers will be rounded off to integers.
#'@param method a character object; defines the method used for the calculation
#'of the forecasting intervals; with \code{"norm"} the intervals are obtained
#'under the assumption of normally distributed innovations; with \code{"boot"}
#'the intervals are obtained via a bootstrap; is set to \code{"norm"} by
#'default.
#'@param alpha a numeric vector of length 1 with \eqn{0 < } \code{alpha}
#'\eqn{ < 1}; the forecasting intervals will be obtained based on the
#'confidence level (\eqn{100}\code{alpha})-percent; is set to
#'\code{alpha = 0.95} by default, i.e., a \eqn{95}-percent confidence level.
#'@param it an integer that represents the total number of iterations, i.e.,
#'the number of simulated series; is set to \code{10000} by default; only
#'necessary, if \code{method = "boot"}; decimal
#'numbers will be rounded off to integers.
#'@param n.start an integer that defines the 'burn-in' number
#'of observations for the simulated ARMA series via bootstrap; is set to
#'\code{1000} by default; only necessary, if \code{method = "boot"};decimal
#'numbers will be rounded off to integers.
#'@param pb a logical value; for \code{pb = TRUE}, a progress bar will be shown
#'in the console, if \code{method = "boot"}.
#'@param cores an integer value >0 that states the number of (logical) cores to
#'use in the bootstrap (or \code{NULL}); the default is the maximum number of
#'available cores
#'(via \code{\link[future:availableCores]{future::availableCores}}); for
#'\code{cores = NULL}, parallel computation is disabled.
#'@param np.fcast a character object; defines the forecasting method used
#'for the nonparametric trend; for \code{np.fcast = "lin"} the trend is
#'is extrapolated linearly based on the last two trend estimates; for
#'\code{np.fcast = "const"}, the last trend estimate is used as a constant
#'estimate for future values; is set to \emph{"lin"} by default.
#'@param export.error a single logical value; if the argument is set to
#'\code{TRUE} and if also \code{method = "boot"}, a list is returned instead
#'of a matrix (\code{FALSE}); the first element of the list is the usual
#'forecasting matrix whereas the second element is a matrix with \code{h}
#'columns, where each column represents the calculated forecasting errors for
#'the respective future time point \eqn{n + 1, n + 2, ..., n + h}; is set to
#'\code{FALSE} by default.
#'@param plot a logical value that controls the graphical output; for
#'\code{plot = TRUE}, the original series with the obtained point forecasts
#'as well as the forecasting intervals will be plotted; for the default
#'\code{plot = FALSE}, no plot will be created.
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
#'This function is part of the \emph{smoots} package and was implemented under
#'version 1.1.0. The point forecasts and forecasting intervals are obtained
#'based on the additive nonparametric regression model
#'\deqn{y_t = m(x_t) + \epsilon_t,}{y_[t] = m(x_[t]) + \epsilon_[t],}
#'where \eqn{y_t}{y_[t]} is the observed time series with equidistant design,
#'\eqn{x_t}{x_[t]} is the rescaled time on the interval \eqn{[0, 1]},
#'\eqn{m(x_t)}{m(x_[t])} is a smooth trend function and
#'\eqn{\epsilon_t}{\epsilon_[t]} are stationary errors with
#'\eqn{E(\epsilon_t) = 0}{E(\epsilon_[t]) = 0} and
#'short-range dependence (see also Beran and Feng, 2002). Thus, we assume
#'\eqn{y_t}{y_[t]} to be a trend-stationary time series. Furthermore, we assume
#'that the rest term \eqn{\epsilon_t}{\epsilon_[t]} follows an ARMA(\eqn{p,q})
#'model
#'\deqn{\epsilon_t = \zeta_t + \beta_1 \epsilon_{t-1} + ... + \beta_p
#'\epsilon_{t-p} + \alpha_1 \zeta_{t-1} + ... +
#'\alpha_q \zeta_{t-q},}{\epsilon_[t] = \zeta_[t] + \beta_[1] \epsilon_[t-1]
#'+ ... + \beta_[p] \epsilon_[t-p] + \alpha_[1] \zeta_[t-1] + ... + \alpha_[q]
#'\zeta_[t-q],}
#'where \eqn{\alpha_j}{\alpha_[j]}, \eqn{j = 1, 2, ..., q}, and
#'\eqn{\beta_i}{\beta_[i]}, \eqn{i = 1, 2, ..., p}, are real numbers and
#'the random variables \eqn{\zeta_t}{\zeta_[t]} are
#'i.i.d. (identically and independently distributed) with
#'zero mean and constant variance.
#'
#'The point forecasts and forecasting intervals for the future periods
#'\eqn{n + 1, n + 2, ..., n + h} will be obtained. With respect to the point
#'forecasts of \eqn{\epsilon_t}{\epsilon_[t]}, i.e.,
#'\eqn{\hat{\epsilon}_{n+k}}{hat[\epsilon]_[n+k]}, where
#'\eqn{k = 1, 2, ..., h},
#'\deqn{\hat{\epsilon}_{n+k} = \sum_{i=1}^{p} \hat{\beta}_i \epsilon_{n+k-i} +
#'\sum_{j=1}^{q} \hat{\alpha}_j \hat{\zeta}_{n+k-j}}{hat[\epsilon]_[n+k] =
#'sum_[i=1]^[p] \{hat[\beta]_[i] \epsilon_[n+k-i]\} +
#'sum_[j=1]^[q] \{hat[\alpha]_[j] hat[\zeta]_[n+k-j]\}}
#'with \eqn{\epsilon_{n+k-i} = \hat{\epsilon}_{n+k-i}}{\epsilon_[n+k-i] =
#'hat[\epsilon]_[n+k-i]} for \eqn{n+k-i > n} and
#'\eqn{\hat{\zeta}_{n+k-j} = E(\zeta_t) = 0}{hat[\zeta]_[n+k-j] = E(\zeta_[t])
#'= 0} for \eqn{n+k-j > n} will be applied. In practice, this procedure will
#'not be applied directly to \eqn{\epsilon_t}{\epsilon_[t]} but to
#'\eqn{y_t - \hat{m}(x_t)}{y_[t] - hat[m](x_[t])}.
#'
#'The point forecasts of the nonparametric trend are simply obtained following
#'the proposal by Fritz et al. (forthcoming) by
#'\deqn{\hat{m}(x_{n+k}) = \hat{m}(x_n) + Dk(\hat{m}(x_n) -
#'\hat{m}(x_{n-1})),}{hat[m](x_[n+k]) = hat[m](x_[n]) + D * k(hat[m](x_[n]) -
#'hat[m](x_[n-1])),}
#'where \eqn{D} is a dummy variable that is either equal to the constant value
#'\eqn{1} or \eqn{0}. Consequently, if \eqn{D = 0},
#'\eqn{\hat{m}(x_{n})}{hat[m](x_[n])}, i.e., the last trend estimate, is
#'used as a constant estimate for the future. However, if \eqn{D = 1}, the
#'trend is extrapolated linearly. The point forecast for the whole component
#'model is then given by
#'\deqn{\hat{y}_{n+k} = \hat{m}(x_{n+k}) + \hat{\epsilon}_{n+k},}{hat[y]_[n+k]
#'= hat[m](x_[n+k]) + hat[\epsilon]_[n+k],}
#'i.e., it is equal to the sum of the point forecasts of the individual
#'components.
#'
#'Equivalently to the point forecasts, the forecasting intervals are the sum
#'of the forecasting intervals of the individual components. To simplify the
#'process, the forecasting error in \eqn{\hat{m}(x_{n+k})}{hat[m](x_[n+k])},
#'which is of order \eqn{O(-2/5)}, is not considered (see Fritz et al.
#'(forthcoming)), i.e., only the forecasting intervals with respect to the
#'rest term \eqn{\epsilon_t}{\epsilon_[t]} will be calculated.
#'
#'If the distribution of the innovations is non-normal or generally not further
#'specified, bootstrapping the forecasting intervals is recommended. If they
#'are however normally distributed or if it is at least assumed that they are,
#'the forecasting errors are also approximately normally distributed with a
#'quickly obtainable variance. For further details on the bootstrapping
#'method, we refer the readers to \code{\link{bootCast}}, whereas more
#'information on the calculation under normality can be found at
#'\code{\link{normCast}}.
#'
#'In order to apply the function, a \code{smoots} object that was generated as
#'the result of a trend estimation process needs to be passed to the argument
#'\code{obj}. The arguments \code{p} and \code{q} represent the orders of the
#'of the ARMA(\eqn{p,q}) model that the error term
#'\eqn{\epsilon_t}{\epsilon_[t]} is assumed to follow. If both arguments are
#'set to \code{NULL}, which is the default setting, orders will be selected
#'according to the Bayesian Information Criterion (BIC) for all possible
#'combinations of \eqn{p,q = 0, 1, ..., 5}. Furthermore, the forecasting
#'horizon can be adjusted by means of the argument \code{h}, so that point
#'forecasts and forecasting intervals will be obtained for all time points
#'\eqn{n + 1, n + 2, ..., n + h}.
#'
#'The function also allows for two calculation approaches for the forecasting
#'intervals. Via the argument \code{method}, intervals
#'can be obtained under the assumption that the ARMA innovations are normally
#'distributed (\code{method = "norm"}). Alternatively, bootstrapped intervals
#'can be obtained for unknown innovation distributions that are clearly
#'non-Gaussian (\code{method = "boot"}).
#'
#'Another argument is \code{alpha}. By passing a value
#'to this argument, the (\eqn{100}\code{alpha})-percent confidence level for
#'the forecasting intervals can be defined. If \code{method = "boot"} is
#'selected, the additional arguments \code{it} and \code{n.start} can be
#'adjusted. More specifically, \code{it} regulates the number of iterations of
#'the bootstrap, whereas \code{n.start} sets the number of 'burn-in'
#'observations in the simulated ARMA processes within the bootstrap that are
#'omitted.
#'
#'Since this bootstrap approach for \code{method = "boot"} generally needs a
#'lot of computation time, especially for
#'series with high numbers of observations and when fitting models with many
#'parameters, parallel computation of the bootstrap iterations is enabled.
#'With \code{cores}, the number of cores can be defined with an integer.
#'Nonetheless, for \code{cores = NULL}, no cluster is created and therefore
#'the parallel computation is disabled. Note that the bootstrapped results are
#'fully reproducible for all cluster sizes. The progress of the bootstrap can
#'be observed in the R console, where a progress bar and the estimated
#'remaining time are displayed for \code{pb = TRUE}.
#'
#'Moreover, the argument \code{np.fcast} allows to set the forecasting method
#'for the nonparametric trend function. As previously discussed, the two
#'options are a linear extrapolation of the trend (\code{np.fcast = "lin"}) and
#'a constant continuation of the last estimated value of the trend
#'(\code{np.fcast = "const"}).
#'
#'The function also implements the option to automatically create a plot of
#'the forecasting results for \code{plot = TRUE}. This includes the feature
#'to pass additional arguments of the standard plot function to
#'\code{modelCast} (see also the section 'Examples').
#'
#'NOTE:
#'
#'Within this function, the \code{\link[stats]{arima}} function of the
#'\code{stats} package with its method \code{"CSS-ML"} is used throughout
#'for the estimation of ARMA models. Furthermore, to increase the performance,
#'C++ code via the \code{\link[Rcpp:Rcpp-package]{Rcpp}} and
#'\code{\link[RcppArmadillo:RcppArmadillo-package]{RcppArmadillo}} packages was
#'implemented. Also, the \code{\link[future:multisession]{future}} and
#'\code{\link[future.apply:future.apply-package]{future.apply}} packages are
#'considered for parallel computation of bootstrap iterations. The progress
#'of the bootstrap is shown via the
#'\code{\link[progressr:progressr-package]{progressr}} package.
#'
#'@md
#'
#'@export
#'
#'@importFrom graphics polygon lines points
#'@importFrom utils tail
#'
#'@return The function returns a \eqn{3} by \eqn{h} matrix with its columns
#'representing the future time points and the point forecasts, the lower
#'bounds of the forecasting intervals and the upper bounds of the
#'forecasting intervals as the rows. If the argument \code{plot} is set to
#'\code{TRUE}, a plot of the forecasting results is created.
#'
#'#'If \code{export.error = TRUE} is selected, a list with the following
#'elements is returned instead.
#'
#'\describe{
#'\item{fcast}{the \eqn{3} by \eqn{h} forecasting matrix with point forecasts
#'and bounds of the forecasting intervals.}
#'\item{error}{an \code{it} by \eqn{h} matrix, where each column
#'represents a future time point \eqn{n + 1, n + 2, ..., n + h}; in each column
#'the respective \code{it} simulated forecasting errors are saved.}
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
#' Feng, Y., Gries, T., Fritz, M., Letmathe, S. and Schulz, D. (2020).
#' Diagnosing the trend and bootstrapping the forecasting intervals using a
#' semiparametric ARMA. Discussion Paper. Paderborn University. Unpublished.
#'
#' Fritz, M., Forstinger, S., Feng, Y., and Gries, T. (forthcoming).
#' Forecasting economic growth processes for developing economies.
#' Unpublished.
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
#'\donttest{
#'X <- log(smoots::gdpUS$GDP)
#'NPest <- smoots::msmooth(X)
#'modelCast(NPest, h = 5, plot = TRUE, xlim = c(261, 295), type = "b",
#'  col = "deepskyblue4", lty = 3, pch = 20, main = "Exemplary title")
#'}
#'

modelCast <- function(obj, p = NULL, q = NULL, h = 1,
  method = c("norm", "boot"), alpha = 0.95, it = 10000, n.start = 1000,
  pb = TRUE, cores = future::availableCores(),
  np.fcast = c("lin", "const"), export.error = FALSE, plot = FALSE,
  ...) {

  if (!inherits(obj, "smoots") || !(attr(obj, "function") %in% c("msmooth",
    "tsmooth", "gsmooth", "knsmooth")) ||
    (attr(obj, "function") == "gsmooth" && obj$v != 0)) {
    stop("The argument 'obj' must be an object of class 'smoots' and the result
         of a trend estimation.")
  }
  if (!(length(method) %in% c(1, 2)) || !all(!is.na(method)) ||
      !is.character(method)) {
    stop("The argument 'method' is defined incorrectly.")
  }
  if (all(method == c("norm", "boot"))) method <- "norm"
  if (length(method) != 1 || !(method %in% c("norm", "boot"))) {
    stop("The argument 'method' is defined incorrectly.")
  }
  if (length(plot) != 1 || is.na(plot) || !(plot %in% c(TRUE, FALSE))) {
    stop("The argument 'plot' must be a single logical value (TRUE or FALSE).")
  }

  trend.fc <- trendCast(object = obj, h = h, np.fcast = np.fcast, plot = FALSE)
  X <- obj$res

  if (is.null(p) && is.null(q)) {
    message("Model selection in progress.")
    p.max <- 5
    q.max <- 5
    BICmat <- unname(critMatrix(X, p.max = p.max, q.max = q.max,
      include.mean = FALSE, criterion = "bic"))
    pq.opt <- which(BICmat == min(BICmat), arr.ind = TRUE) - 1
    p <- pq.opt[[1]]
    q <- pq.opt[[2]]
    message("Orders p=", p, " and q=", q, " were selected.")
  }
  if (is.null(p)) p <- 0
  if (is.null(q)) q <- 0
  if (method == "norm") {
    FI <- normCast(X = X, p, q, include.mean = FALSE,
      h = h, alpha = alpha, plot = FALSE)
  } else if (method == "boot" && export.error == FALSE) {
    FI <- bootCast(X = X, p, q, include.mean = FALSE,
      n.start = n.start, h = h, it = it, alpha = alpha, pb = pb,
      cores = cores, export.error = FALSE, plot = FALSE)
  } else if (method == "boot" && export.error == TRUE) {
    bootCalculation <- bootCast(X = X, p, q, include.mean = FALSE,
      n.start = n.start, h = h, it = it, alpha = alpha, pb = pb,
      cores = cores, export.error = TRUE, plot = FALSE)
    FI <- bootCalculation[["fcast"]]
  }

  for (i in 1:3) {
    FI[i, ] <- FI[i, ] + trend.fc
  }

  if (plot == TRUE) {
    n <- length(X)
    dots <- list(...)
    if (is.null(dots[["x"]])) dots[["x"]] <- 1:n
    dots[["y"]] <- X <- obj$orig
    if (is.null(dots[["xlab"]])) dots[["xlab"]] <- "Time"
    if (is.null(dots[["ylab"]])) dots[["ylab"]] <- "Series"
    if (is.null(dots[["main"]])) dots[["main"]] <- "Forecasting results"
    if (is.null(dots[["type"]])) dots[["type"]] <- "l"
    xStep <- dots$x[2] - dots$x[1]
    x.fc <- seq(from = dots$x[n] + xStep, to = dots$x[n] + h * xStep,
      by = xStep)
    if (is.null(dots[["xlim"]])) {
      xlim1 <- tail(dots$x, 6 * h)[1]
      dots[["xlim"]] <- c(xlim1, dots$x[n] + h * xStep)
    }
    if (is.null(dots[["ylim"]])) {
      x.all <- c(dots$x, x.fc)
      low.bound <- c(obj[["orig"]], FI[2, ])
      up.bound <- c(obj[["orig"]], FI[3, ])
      x.low <- sum(x.all < dots$xlim[1]) + 1
      x.up <- sum(x.all <= dots$xlim[2])
      dots[["ylim"]] <- c(min(low.bound[x.low:x.up]),
        max(up.bound[x.low:x.up]))
    }
    do.call(what = graphics::plot, args = dots)
    if (h >= 2) {
      polygon(c(x.fc, rev(x.fc)),
        c(FI[2, ], rev(FI[3, ])), col = col.alpha("gray", alpha = 0.8),
        border = NA)
      lines(x.fc, FI[1, ], col = "red", lty = 1)
    } else {
      lines(c(x.fc, x.fc), c(FI[2, ], FI[3, ]),
        col = col.alpha("gray", alpha = 0.8), lty = 1)
      points(x.fc, FI[1, ], col = "red", pch = 20)
    }
  }

  if (method == "boot" && export.error == TRUE) {
    FI <- list(fcast = FI, error = bootCalculation[["error"]])
  }

  FI
}
