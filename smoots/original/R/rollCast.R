#' Backtesting Semi-ARMA Models with Rolling Forecasts
#'
#' A simple backtest of Semi-ARMA models via rolling forecasts can be
#' implemented.
#'
#'@param y a numeric vector that represents the equidistant time series assumed
#'to follow a Semi-ARMA model; must be ordered from past to present.
#'@param p an integer value \eqn{\geq 0}{\ge 0} that defines the AR order
#'\eqn{p} of the underlying ARMA(\eqn{p,q}) model within \code{X}; is set to
#'\code{NULL} by default; if no value is passed to \code{p} but one is passed
#'to \code{q}, \code{p} is set to \code{0}; if both \code{p} and \code{q} are
#'\code{NULL}, optimal orders following the BIC for
#'\eqn{0 \leq p,q \leq 5}{0 \le p,q \le 5} are chosen; is set to \code{NULL} by
#'default; decimal numbers will be rounded off to integers.
#'@param q an integer value \eqn{\geq 0}{\ge 0} that defines the MA order
#'\eqn{q} of the underlying ARMA(\eqn{p,q}) model within \code{X}; is set to
#'\code{NULL} by default; if no value is passed to \code{q} but one is passed
#'to \code{p}, \code{q} is set to \code{0}; if both \code{p} and \code{q} are
#'\code{NULL}, optimal orders following the BIC for
#'\eqn{0 \leq p,q \leq 5}{0 \le p,q \le 5} are chosen; is set to \code{NULL} by
#'default; decimal numbers will be rounded off to integers.
#'@param K a single, positive integer value that defines the number of
#'out-of-sample observations; the last \code{K} observations in \code{y} are
#'treated as the out-of-sample observations, whereas the rest of the
#'observations in \code{y} are the in-sample values.
#'@param method a character object; defines the method used for the calculation
#'of the forecasting intervals; with \code{"norm"} the intervals are obtained
#'under the assumption of normally distributed innovations; with \code{"boot"}
#'the intervals are obtained via a bootstrap; is set to \code{"norm"} by
#'default.
#'@param alpha a numeric vector of length 1 with \eqn{0 < } \code{alpha}
#'\eqn{ < 1}; the forecasting intervals will be obtained based on the
#'confidence level (\eqn{100}\code{alpha})-percent; is set to
#'\code{alpha = 0.95} by default, i.e., a \eqn{95}-percent confidence level.
#'@param np.fcast a character object; defines the forecasting method used
#'for the nonparametric trend; for \code{np.fcast = "lin"} the trend is
#'is extrapolated linearly based on the last two trend estimates; for
#'\code{np.fcast = "const"}, the last trend estimate is used as a constant
#'estimate for future values; is set to \emph{"lin"} by default.
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
#'@param argsSmoots a list that contains arguments that will be passed to
#'\code{\link{msmooth}} for the estimation of the nonparametric trend
#'function; by default, the default values of \code{msmooth} are used.
#'@param plot a logical value that controls the graphical output; for the
#'default (\code{plot = TRUE}), the original series with the obtained point
#'forecasts as well as the forecasting intervals will be plotted; for
#'\code{plot = FALSE}, no plot will be created.
#'@param argsPlot a list; additional arguments for the standard plot function,
#'e.g., \code{xlim}, \code{type}, ..., can be passed to it; arguments with
#'respect to plotted graphs, e.g., the argument \code{col}, only affect the
#'original series \code{y}; please note that in accordance with the argument
#'\code{x} (lower case) of the standard plot function, an additional numeric
#'vector with time points can be implemented via the argument \code{x} (lower
#'case).
#'
#'@export
#'
#'@importFrom stats arima qnorm quantile
#'@importFrom utils head tail
#'@importFrom graphics lines points polygon
#'
#'@details
#'Define that an observed, equidistant time series \eqn{y_t}{y_[t]}, with
#'\eqn{t = 1, 2, ..., n}, follows
#'\deqn{y_t = m(x_t) + \epsilon_t,}{y_[t] = m(x_[t]) + \epsilon_[t],}
#'where \eqn{x_t = t/n}{x_[t] = t/n} is the rescaled time on the closed
#'interval \eqn{[0,1]} and \eqn{m(x_t)}{m(x_[t])} is a nonparametric and
#'deterministic trend function (see Beran and Feng, 2002, and Feng, Gries and
#'Fritz, 2020).
#'\eqn{\epsilon_t}{\epsilon_[t]}, on the other hand, is a stationary process
#'with \eqn{E(\epsilon_t) = 0}{E(\epsilon_[t]) = 0} and short-range dependence.
#'For the purpose of this function, \eqn{\epsilon_t}{\epsilon_[t]} is assumed
#'to follow an autoregressive-moving-average (ARMA) model with
#'\deqn{\epsilon_t = \zeta_t + \beta_1 \epsilon_{t-1} + ... + \beta_p
#'\epsilon_{t-p} + \alpha_1 \zeta_{t-1} + ... +
#'\alpha_q \zeta_{t-q}.}{\epsilon_[t] = \zeta_[t] + \beta_[1] \epsilon_[t-1] +
#'... + \beta_[p] \epsilon_[t-p] + \alpha_[1] \zeta_[t-1] + ... +
#'\alpha_[q] \zeta_[t-q].}
#'Here, the random variables \eqn{\zeta_t}{\zeta_[t]} are identically and
#'independently distributed (i.i.d.) with zero-mean and a constant variance
#'and the coefficients \eqn{\alpha_j}{\alpha_[j]} and \eqn{\beta_i}{\beta_[i]},
#'\eqn{i = 1, 2, ..., p} and \eqn{j = 1, 2, ..., q}, are real numbers.
#'The combination of both previous formulas will be called a semiparametric
#'ARMA (Semi-ARMA) model.
#'
#'An explicit forecasting method of Semi-ARMA models is described in
#'\code{\link{modelCast}}. To backtest a selected model, a slightly adjusted
#'procedure is used. The data is divided into in-sample and an
#'out-of-sample values (usually the last \eqn{K = 5} observations in the data
#'are reserved for the out-of-sample observations). A model is fitted to the
#'in-sample data, whereas one-step rolling point forecasts and forecasting
#'intervals are obtained for the out-of-sample time points. The proposed
#'forecasts of the trend are either a linear or a constant extrapolation of
#'the trend with negligible forecasting intervals, whereas the point forecasts
#'of the stationary rest term are obtained via the selected ARMA(\eqn{p,q})
#'model (see Fritz et al., 2020). The corresponding forecasting intervals
#'are calculated under the assumption that the innovations
#'\eqn{\zeta_t}{\zeta_[t]} are either normally distributed (see e.g. pp.
#'93-94 in Brockwell and Davis, 2016) or via a forward bootstrap (see Lu and
#'Wang, 2020). For a one-step forecast for time point \eqn{t}, all observations
#'until time point \eqn{t-1} are assumed to be known.
#'
#'The function calculates three important values for backtesting: the number
#'of breaches, i.e. the number of true observations that lie outside of the
#'forecasting intervals, the mean absolute scaled error (MASE, see Hyndman
#'and Koehler, 2006) and the root mean squared scaled error (RMSSE, see
#'Hyndman and Koehler, 2006) are obtained. For the MASE, a value \eqn{< 1}
#'indicates a better average forecasting potential than a naive forecasting
#'approach.
#'Furthermore, it is independent from the scale of the data and can thus be
#'used to compare forecasts of different datasets. Closely related is the
#'RMSSE, however here, the mean of the squared forecasting errors is computed
#'and scaled by the mean of the squared naive forecasting approach. Then the
#'root of that value is the RMSSE. Due to the close relation, the
#'interpretation of the RMSSE is similarly but not identically to the
#'interpretation of the MASE. Of course, a value close to zero is preferred
#'in both cases.
#'
#'To make use of the function, a numeric vector with the values of a time
#'series that is assumed to follow a Semi-ARMA model needs to be passed to
#'the argument \code{y}. Moreover, the arguments \code{p} and \code{q}
#'represent the AR and MA orders, respectively, of the underlying ARMA
#'process in the parametric part of the model. If both values are set to
#'\code{NULL}, an optimal order in accordance with the Bayesian Information
#'Criterion (BIC) will be selected. If only one of the values is \code{NULL},
#'it will be changed to zero instead. \code{K} defines the number of the
#'out-of-sample observations; these will be cut off the end of \code{y}, while
#'the remaining observations are treated as the in-sample observations. For the
#'\eqn{K} out-of-sample time points, rolling forecasts will be obtained.
#'\code{method} describes the method to use for the computation of the
#'prediction intervals. Under the normality assumption for the innovations
#'\eqn{\zeta_t}{\zeta_[t]}, intervals can be obtained via
#'\emph{method} = "norm". However, if the assumption does not hold, a
#'bootstrap can be implemented as well (\emph{method = "boot"}). Both
#'approaches are explained in more detail in \code{\link{normCast}} and
#'\code{\link{bootCast}}, respectively. With \code{alpha}, the confidence
#'level of the forecasting intervals can be adjusted, as the
#'(\eqn{100}\code{alpha})-percent forecasting intervals will be computed. By
#'means of the argument \code{np.fcast}, the forecasting method for the
#'nonparametric trend function can be defined. Selectable are a linear
#'(\code{np.fcast = "lin"}) and a constant (\code{np.fcast = "const"})
#'extrapolation. For more information on these methods, we refer the reader to
#'\code{\link{trendCast}}.
#'
#'\code{it}, \code{n.start}, \code{pb} and \code{cores} are only
#'relevant for \code{method = "boot"}. With \code{it} the total number of
#'bootstrap iterations is defined, whereas \code{n.start} regulates, how
#'many 'burn-in' observations are generated for each simulated ARMA process
#'in the bootstrap. Since a bootstrap may take a longer computation time,
#'with the argument \code{cores} the number of cores for parallel computation
#'of the bootstrap iterations can be defined. Nonetheless, for
#'\code{cores = NULL}, no cluster is created and therefore
#'the parallel computation is disabled. Note that the bootstrapped results are
#'fully reproducible for all cluster sizes. Moreover, for \code{pb = TRUE},
#'the progress of the bootstrap approach can be observed in the R console via
#'a progress bar. Additional information on these four function arguments can
#'be found in \code{\link{bootCast}}.
#'
#'The argument \code{argsSmoots} is a list. In this list, different arguments
#'of the function \code{\link{msmooth}} can be implemented to adjust the
#'estimation of the nonparametric part of the complete model. The arguments
#'of the smoothing function are described in \code{\link{msmooth}}.
#'
#'\code{rollCast} allows for a quick plot of the results. If the logical
#'argument \code{plot} is set to \code{TRUE}, a graphic with default
#'settings is created. Nevertheless, users are allowed to implement further
#'arguments of the standard plot function in the list \code{argsPlot}. For
#'example, the limits of the plot can be adjusted by \code{xlim} and
#'\code{ylim}. Furthermore, an argument \code{x} can be included in
#'\code{argsPlot} with the actual equidistant time points of the whole series
#'(in-sample & out-of-sample observations). Otherwise, simply \code{1:n} is
#'used as the in-sample time points by default.
#'
#'NOTE:
#'
#'Within this function, the \code{\link[stats]{arima}} function of the
#'\code{stats} package with its method \code{"CSS-ML"} is used throughout for
#'the estimation of ARMA models. Furthermore, to increase the performance,
#'C++ code via the \code{\link[Rcpp:Rcpp-package]{Rcpp}} and
#'\code{\link[RcppArmadillo:RcppArmadillo-package]{RcppArmadillo}} packages
#'was implemented. Also, the \code{\link[future:multisession]{future}} and
#'\code{\link[future.apply:future.apply-package]{future.apply}} packages are
#'considered for parallel computation of bootstrap iterations. The progress
#'of the bootstrap is shown via the
#'\code{\link[progressr:progressr-package]{progressr}} package.
#'
#'@md
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
#'@return
#'A list with different elements is returned. The elements are as follows.
#'\describe{
#'\item{alpha}{a single numeric value; it describes, what confidence level
#'(\eqn{100}\code{alpha})-percent has been considered for the forecasting
#'intervals.}
#'\item{breach}{a logical vector that states whether the \eqn{K} true
#'out-of-sample observations lie outside of the forecasting intervals,
#'respectively; a breach is denoted by \code{TRUE}.}
#'\item{breach.val}{a numeric vector that contains the margin of the breaches
#'(in absolute terms) for the \eqn{K} out-of-sample time points; if a breach
#'did not occur, the respective element is set to zero.}
#'\item{error}{a numeric vector that contains the simulated empirical
#'values of the forecasting error for \code{method = "boot"}; otherwise,
#'it is set to \code{NULL}.}
#'\item{fcast.rest}{a numeric vector that contains the \eqn{K} point forecasts
#'of the parametric part of the model.}
#'\item{fcast.roll}{a numeric matrix that contains the \eqn{K} rolling point
#'forecasts as well as the values of the bounds of the respective forecasting
#'intervals for the complete model;
#'the first row contains the point forecasts, the lower bounds of the
#'forecasting intervals are in the second row and the upper bounds
#'can be found in the third row.}
#'\item{fcast.trend}{a numeric vector that contains the \eqn{K} obtained trend
#'forecasts.}
#'\item{K}{a positive integer; states the number of out-of-sample observations
#'as well as the number of forecasts for the out-of-sample time points.}
#'\item{MASE}{the obtained value of the mean average scaled error for the
#'selected model.}
#'\item{method}{a character object that states, whether the forecasting
#'intervals were obtained via a bootstrap (\code{method = "boot"}) or under
#'the normality assumption for the innovations (\code{method = "norm"}).}
#'\item{model.nonpar}{the output (usually a list) of the nonparametric
#'trend estimation via \code{\link{msmooth}}.}
#'\item{model.par}{the output (usually a list) of the parametric ARMA
#'estimation of the detrended series via \code{\link[stats]{arima}}.}
#'\item{n}{the number of observations (in-sample & out-of-sample
#'observations).}
#'\item{n.in}{the number of in-sample observations (\code{n - n.out}).}
#'\item{n.out}{the number of out-of-sample observations (equals \code{K}).}
#'\item{np.fcast}{a character object that states the applied forecasting
#'method for the nonparametric trend function; either a linear (
#'\code{np.fcast = "lin"}) or a constant \code{np.fcast = "const"} are
#'possible.}
#'\item{quants}{a numeric vector of length 2 with the
#'\eqn{[100(1 -} \code{alpha}\eqn{)/2]}-percent and
#'\{\eqn{100}\eqn{[1 - (1 -} \code{alpha}\eqn{)/2]}\}-percent quantiles of
#'the forecasting error distribution.}
#'\item{RMSSE}{the obtained value of the root mean squared scaled error for
#'the selected model.}
#'\item{y}{a numeric vector that contains all true observations (in-sample &
#'out-of-sample observations).}
#'\item{y.in}{a numeric vector that contains all in-sample observations.}
#'\item{y.out}{a numeric vector that contains the \eqn{K} out-of-sample
#'observations.}
#'}
#'
#'@references
#'Beran, J., and Feng, Y. (2002). Local polynomial fitting with long-memory,
#'short-memory and antipersistent errors. Annals of the Institute of
#'Statistical Mathematics, 54, 291-311.
#'
#'Brockwell, P. J., and Davis, R. A. (2016). Introduction to time series
#'and forecasting, 3rd edition. Springer.
#'
#'Fritz, M., Forstinger, S., Feng, Y., and Gries, T. (2020). Forecasting
#'economic growth processes for developing economies. Unpublished.
#'
#'Feng, Y., Gries, T. and Fritz, M. (2020). Data-driven
#'local polynomial for the trend and its derivatives in economic time
#'series. Journal of Nonparametric Statistics, 32:2, 510-533.
#'
#'Hyndman, R. J., and Koehler, A. B. (2006). Another look at measures of
#'forecast accuracy. International Journal of Forecasting, 22:4, 679-688.
#'
#'Lu, X., and Wang, L. (2020). Bootstrap prediction interval for ARMA models
#'with unknown orders. REVSTAT–Statistical Journal, 18:3, 375-396.
#'
#'@examples
#'lgdp <- log(smoots::gdpUS$GDP)
#'time <- seq(from = 1947.25, to = 2019.5, by = 0.25)
#'backtest <- rollCast(lgdp, K = 5,
#'  argsPlot = list(x = time, xlim = c(2012, 2019.5), col = "forestgreen",
#'  type = "b", pch = 20, lty = 2, main = "Example"))
#'backtest
#'

rollCast <- function(y, p = NULL, q = NULL, K = 5, method = c("norm", "boot"),
  alpha = 0.95, np.fcast = c("lin", "const"), it = 10000, n.start = 1000,
  pb = TRUE, cores = future::availableCores(),
  argsSmoots = list(), plot = TRUE, argsPlot = list()) {

  if (length(y) <= 1 || !all(!is.na(y)) || !is.numeric(y)) {
    stop("A numeric vector without NAs must be passed to the argument 'y'.")
  }

  if (length(K) != 1 || is.na(K) || !is.numeric(K) || K == 0) {
    stop("A single integer value > 0 must be passed to the argument 'K'.")
  }
  K <- floor(K)

  if (!(length(method) %in% c(1, 2)) || !all(!is.na(method)) ||
      !is.character(method)) {
    stop("The argument 'method' is defined incorrectly.")
  }
  if (all(method == c("norm", "boot"))) method <- "norm"
  if (length(method) != 1 || !(method %in% c("norm", "boot"))) {
    stop("The argument 'method' is defined incorrectly.")
  }

  if (!(length(np.fcast) %in% c(1, 2)) || !all(!is.na(np.fcast)) ||
      !is.character(np.fcast)) {
    stop("Argument 'np.fcast' is defined incorrectly.")
  }
  if (all(np.fcast == c("lin", "const"))) np.fcast <- "lin"
  if (length(np.fcast) != 1 || !(np.fcast %in% c("lin", "const"))) {
    stop("Argument 'np.fcast' is defined incorrectly.")
  }

  if (length(plot) != 1 || is.na(plot) || !is.logical(plot)) {
    stop("A single logical value muste be passed to the argument 'plot'.")
  }

  n <- length(y)
  n.out <- K
  n.in <- n - n.out
  y.in <- y[1:n.in]
  y.out <- y[(n.in + 1):(n)]

  argsSmoots[["y"]] <- y.in
  np.est <- suppressMessages(do.call(what = msmooth,
    args = argsSmoots))

  if (is.null(p) && is.null(q)) {
    message("Model selection in progress.")
    p.max <- 5
    q.max <- 5
    BICmat <- unname(critMatrix(np.est[["res"]], p.max = p.max, q.max = q.max,
       include.mean = FALSE, criterion = "bic"))
    pq.opt <- which(BICmat == min(BICmat), arr.ind = TRUE) - 1
    p <- pq.opt[[1]]
    q <- pq.opt[[2]]
    message("Orders p=", p, " and q=", q, " were selected.")
  }
  if (is.null(p)) p <- 0
  if (is.null(q)) q <- 0
  p <- floor(p)
  q <- floor(q)
  if (np.fcast == "lin") {
    trendf <- np.est[["ye"]][[n.in]] - np.est[["ye"]][[n.in - 1]]
  } else {
    trendf <- 0
  }
  fcast.trend <- np.est[["ye"]][[n.in]] + (1:K) * trendf
  arma.est <- suppressWarnings(arima(np.est[["res"]],
    order = c(p, 0, q), include.mean = FALSE))
  ma <- ar <- 0
  if (p > 0) ar <- unname(arma.est[["coef"]][1:p])
  if (q > 0) ma <- unname(arma.est[["coef"]][(p + 1):(p + q)])
  p.star <- length(ar)
  q.star <- length(ma)

  l <- max(p.star, q.star)
  y.out.arma <- y.out - fcast.trend
  x <- c(tail(np.est[["res"]], n = l), y.out.arma)
  u <- c(tail(arma.est[["residuals"]], n = l), rep(NA, times = K))
  arma.fcast <- rep(NA, times = K)

  for (i in seq_along(arma.fcast)) {
    arma.fcast[i] <- ar %*% x[(i + l - 1):(i + l - p.star)] +
      ma %*% u[(i + l - 1):(i + l - q.star)]
    u[i + l] <- y.out.arma[i] - arma.fcast[i]
  }
  alpha.s <- 1 - alpha
  if (method == "norm") {
    sig <- sqrt(arma.est[["sigma2"]])
    quants <- qnorm(c(alpha.s / 2, 1 - alpha.s / 2), mean = 0, sd = sig)
    names(quants) <- paste0(c(100 * alpha.s / 2, 100 * (1 - alpha.s / 2)), "%")
    fcastErr <- NULL
  } else {
    fcastModel <- bootCast(X = np.est[["res"]], p = p, q = q, h = 1,
      alpha = alpha, it = it, n.start = n.start, pb = pb,
      cores = cores, export.error = TRUE, plot = FALSE)
    fcastErr <- c(fcastModel[["error"]])
    quants <- quantile(fcastErr, probs = c(alpha.s / 2, 1 - alpha.s / 2))
  }

  fcast.compl <- fcast.trend + arma.fcast
  fc.low <- fcast.compl + quants[[1]]
  fc.up <- fcast.compl + quants[[2]]
  fcast <- rbind(fcast.compl, fc.low, fc.up)
  colnames(fcast) <- paste0("k=", 1:K)
  rownames(fcast) <- c("fcast", paste0(c(100 * alpha.s / 2,
    100 * (1 - alpha.s / 2)), "%"))
  breach.up <- y.out > fc.up
  breach.low <- y.out < fc.low
  breach.val <- rep(0, times = K)
  breach.val[breach.up] <- y.out[breach.up] - fc.up[breach.up]
  breach.val[breach.low] <- y.out[breach.low] - fc.low[breach.low]
  breach <- breach.val != 0
  MASE <- mean(abs(y.out - fcast.compl)) / (mean(abs(diff(y.in))))
  RMSSE <- sqrt(mean((y.out - fcast.compl)^2) /
    (mean(diff(y.in)^2)))

  out <- list(fcast.roll = fcast, y.out = y.out, model.nonpar = np.est,
    model.par = arma.est, alpha = alpha, K = K, fcast.trend = fcast.trend,
    fcast.rest = arma.fcast, quants = quants, breach = breach,
    breach.val = breach.val, n = n, n.in = n.in, n.out = n.out,
    method = method, np.fcast = np.fcast, it = it, MASE = MASE, RMSSE = RMSSE,
    error = fcastErr, y = y, y.in = y.in)

  class(out) <- "smoots"
  attr(out, "function") <- "rollCast"

  if (plot == TRUE) {
    if (is.null(argsPlot[["main"]])) {
      argsPlot[["main"]] <- "Backtesting results"
    }
    if (is.null(argsPlot[["xlab"]])) argsPlot[["xlab"]] <- "Time"
    if (is.null(argsPlot[["ylab"]])) argsPlot[["ylab"]] <- "Series"
    if (is.null(argsPlot[["type"]])) argsPlot[["type"]] <- "l"
    if (is.null(argsPlot[["x"]])) argsPlot[["x"]] <- 1:n
    argsPlot[["y"]] <- y.in
    x.out <- tail(argsPlot[["x"]], n.out)
    x.in <- head(argsPlot[["x"]], n.in)
    argsPlot[["x"]] <- x.in
    if (is.null(argsPlot[["xlim"]])) {
      xlim1 <- tail(x.in, 6 * K)[1]
      argsPlot[["xlim"]] <- c(xlim1, tail(x.out, 1))
    }

    if (is.null(argsPlot[["ylim"]])) {
      x.all <- c(x.in, x.out)
      x.low <- sum(x.all < argsPlot[["xlim"]][[1]]) + 1
      x.up <- sum(x.all <= argsPlot[["xlim"]][[2]])
      fc.bound.low <- c(y.in, fc.low)
      fc.bound.up <- c(y.in, fc.up)
      argsPlot[["ylim"]] <- c(min(y[x.low:x.up], fc.bound.low[x.low:x.up]),
        max(y[x.low:x.up], fc.bound.up[x.low:x.up]))
    }

    do.call(what = graphics::plot, args = argsPlot)
    if (K > 1) {
      polygon(c(x.out, rev(x.out)), c(fc.low, rev(fc.up)), border = NA,
        col = col.alpha("gray", 0.8))
      argsPlot[["x"]] <- x.out
      argsPlot[["y"]] <- y.out
      do.call(what = lines, args = argsPlot)
      lines(x.out, fcast.compl, col = "red")
    } else {
      points(c(x.out, x.out), c(fc.up, fc.low), type = "l",
        col = col.alpha("gray", 0.8))
      argsPlot[["type"]] <- "p"
      if (is.null(argsPlot[["pch"]])) argsPlot[["pch"]] <- 20
      argsPlot[["x"]] <- x.out
      argsPlot[["y"]] <- y.out
      do.call(what = points, args = argsPlot)
      points(x.out, fcast.compl, col = "red", pch = 20)
    }
    if (sum(breach >= 1)) {
      points(x.out[breach == TRUE], y.out[breach == TRUE],
        col = "orange")
    }
  }

  out
}
