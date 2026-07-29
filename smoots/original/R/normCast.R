#' Forecasting Function for ARMA Models under Normally Distributed Innovations
#'
#' Point forecasts and the respective forecasting intervals for autoregressive-
#' moving-average (ARMA) models can be calculated, the latter under the
#' assumption of normally distributed innovations, by means of this function.
#'
#'@param X a numeric vector that contains the time series that is assumed to
#'follow an ARMA model ordered from past to present.
#'@param p an integer value \eqn{>= 0} that defines the AR order \eqn{p} of the
#'underlying ARMA(\eqn{p,q}) model within \code{X}; is set to \code{NULL} by
#'default; if no value is passed to \code{p} but one is passed to \code{q},
#'\code{p} is set to \code{0}; if both \code{p} and \code{q} are \code{NULL},
#'optimal orders following the BIC for \eqn{0 \leq p,q \leq 5}{0 \le p,q \le 5}
#'are chosen; is set to \code{NULL} by default; decimal numbers will be rounded
#'off to integers.
#'@param q an integer value \eqn{>= 0} that defines the MA order \eqn{q} of the
#'underlying ARMA(\eqn{p,q}) model within \code{X}; is set to \code{NULL} by
#'default; if no value is passed to \code{q} but one is passed to \code{p},
#'\code{q} is set to \code{0}; if both \code{p} and \code{q} are \code{NULL},
#'optimal orders following the BIC for \eqn{0 \leq p,q \leq 5}{0 \le p,q \le 5}
#'are chosen; is set to \code{NULL} by default; decimal numbers will be rounded
#'off to integers.
#'@param include.mean a logical value; if set to \code{TRUE}, the mean of the
#'series is also estimated; if set to \code{FALSE},
#'\eqn{E(X_t) = 0}{E(X_[t]) = 0} is assumed; is set to \code{FALSE} by default.
#'@param h an integer that represents the forecasting horizon; if \eqn{n} is
#'the number of observations, point forecasts and forecasting intervals will be
#'obtained for the time points \eqn{n + 1} to \eqn{n + h}; is set to \code{1}
#'by default; decimal numbers will be rounded off to integers.
#'@param alpha a numeric vector of length 1 with \eqn{0 <} \code{alpha}
#'\eqn{< 1}; the forecasting intervals will be obtained based on the confidence
#'level (\eqn{100}\code{alpha})-percent; is set to \emph{alpha = 0.95} by
#'default, i.e., a \eqn{95}-percent confidence level.
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
#'\code{length(x) == length(X)} should be \code{TRUE}, as future time
#'points will be calculated automatically.
#'
#'@details
#'This function is part of the \code{smoots} package and was implemented under
#'version 1.1.0. For a given time series \eqn{X_[t]}, \eqn{t = 1, 2, ..., n},
#'the point forecasts and the respective forecasting intervals will be
#'calculated.
#'It is assumed that the series follows an ARMA(\eqn{p,q}) model
#'\deqn{X_t - \mu = \epsilon_t + \beta_1 (X_{t-1} - \mu) + ... + \beta_p
#'(X_{t-p} - \mu) + \alpha_1 \epsilon_{t-1} + ... + \alpha_{q}
#'\epsilon_{t-q},}{X_[t] - \mu = \epsilon_[t] + \beta_[1] (X_[t-1] - \mu) + ...
#'+ \beta_[p] (X_[t-p] - \mu) + \alpha_[1] \epsilon_[t-1] + ... + \alpha_[q]
#'\epsilon_[t-q],}
#'where \eqn{\alpha_j}{\alpha_[j]} and \eqn{\beta_i}{\beta_[i]} are real
#'numbers (for \eqn{i = 1, 2, .., p} and \eqn{j = 1, 2, ..., q}) and
#'\eqn{\epsilon_t}{\epsilon_[t]} are i.i.d. (identically and independently
#'distributed) random variables with zero mean and constant variance.
#'\eqn{\mu} is equal to \eqn{E(X_t)}{E(X_[t])}.
#'
#'The point forecasts and forecasting intervals for the future periods
#'\eqn{n + 1, n + 2, ..., n + h} will be obtained. With respect to the point
#'forecasts \eqn{\hat{X}_{n + k}}{hat[X]_[n+k]}, where \eqn{k = 1, 2, ..., h},
#'\deqn{\hat{X}_{n + k} = \hat{\mu} + \sum_{i = 1}^{p} \hat{\beta}_{i}
#'(X_{n + k - i} - \hat{\mu}) + \sum_{j = 1}^{q} \hat{\alpha}_{j}
#'\hat{\epsilon}_{n + k - j}}{hat[X]_[n+k] = hat[\mu] +
#'sum_[i=1]^[p] \{hat[\beta]_[i](X_[n+k-i] - hat[\mu])\} +
#'sum_[j=1]^[q] \{hat[\alpha]_[j]hat[\epsilon]_[n+k-j]\}}
#'with \eqn{X_{n+k-i} = \hat{X}_{n+k-i}}{X_[n+k-i] = hat[X]_[n+k-i]} for
#'\eqn{n+k-i > n} and
#'\eqn{\hat{\epsilon}_{n+k-j} = E(\epsilon_t) = 0}{hat[\epsilon]_[n+k-j] =
#'E(\epsilon_[t]) = 0} for \eqn{n+k-j > n} will be applied.
#'
#'The forecasting intervals on the other hand are obtained under the assumption
#'of normally distributed innovations. Let \eqn{q(c)} be the \eqn{100c}-percent
#'quantile of the standard normal distribution. The \eqn{100a}-percent
#'forecasting interval at a point \eqn{n + k}, where \eqn{k = 1, 2, ..., h},
#'is given by
#'\deqn{[\hat{X}_{n+k} - q(a_r)s_k, \hat{X}_{n+k} + q(a_r)s_k]}{[hat[X]_[n+k] -
#'q(a_[r])s_[k], hat[X]_[n+k] + q(a_[r])s_[k]]}
#'with \eqn{s_k}{s_[k]} being the standard deviation of the forecasting error
#'at the time point \eqn{n + k} and with
#'\eqn{a_r = 1 - (1 - a)/2}{a_[r] = 1 - (1 - a)/2}. For ARMA models with normal
#'innovations, the variance of the forecasting error can be derived from the
#'MA(\eqn{\infty}{infinity}) representation of  the model. It is
#'\deqn{\sigma_{\epsilon}^{2} \sum_{i=0}^{k - 1} d_{i}^{2},}{\sigma_[\epsilon]
#'^[2] * sum_[i=0]^[k-1] \{d_[i]^[2]\},}
#'where \eqn{d_i}{d_[i]} are the coefficients of the MA(\eqn{\infty}{infinity})
#'representation and \eqn{\sigma_{\epsilon}^{2}}{\sigma_[\epsilon]^[2]} is the
#'innovation variance.
#'
#'The function \code{normCast} allows for different adjustments to
#'the forecasting progress. At first, a vector with the values of the observed
#'time series ordered from past to present has to be passed to the argument
#'\code{X}. Orders \eqn{p} and \eqn{q} of the underlying ARMA process can be
#'defined via the arguments \code{p} and \code{q}. If only one of these orders
#'is inserted by the user, the other order is automatically set to \code{0}. If
#'none of these arguments are defined, the function will choose orders based on
#'the Bayesian Information Criterion (BIC) for \eqn{0 \leq p,q \leq 5}{0 \le
#'p,q \le 5}. Via the logical argument \code{include.mean} the user can decide,
#'whether to consider the mean of the series within the estimation process.
#'Furthermore, the argument \code{h} allows for the definition of the maximum
#'future time point \eqn{n + h}. Point forecasts and forecasting intervals will
#'be returned for the time points \eqn{n + 1} to \eqn{n + h}. Another argument
#'is \code{alpha}, which is the equivalent of the confidence level considered
#'within the calculation of the forecasting intervals, i.e., the quantiles
#'\eqn{(1 -} \code{alpha}\eqn{)/2} and
#'\eqn{1 - (1 -} \code{alpha}\eqn{)/2} of the
#'forecasting intervals will be obtained.
#'
#'For simplicity, the function also incorporates the possibility to directly
#'create a plot of the output, if the argument \code{plot} is set to
#'\code{TRUE}. By the additional and optional arguments \code{...}, further
#'arguments of the standard plot function can be implemented to shape the
#'returned plot.
#'
#'NOTE:
#'Within this function, the \code{\link[stats]{arima}} function of the
#'\code{stats} package with its method \code{"CSS-ML"} is used throughout
#'for the estimation of ARMA models.
#'
#'@return The function returns a \eqn{3} by \eqn{h} matrix with its columns
#'representing the future time points and the point forecasts, the lower
#'bounds of the forecasting intervals and the upper bounds of the
#'forecasting intervals as the rows. If the argument \code{plot} is set to
#'\code{TRUE}, a plot of the forecasting results is created.
#'
#'@export
#'
#'@importFrom stats arima qnorm
#'@importFrom utils tail
#'@importFrom graphics polygon lines points
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
#'\item Dominik Schulz (Research Assistant) (Department of Economics, Paderborn
#'University), \cr
#'Package Creator and Maintainer
#'}
#'
#'@examples
#'### Example 1: Simulated ARMA process ###
#'
#'# Simulation of the underlying process
#'n <- 2000
#'n.start = 1000
#'set.seed(21)
#'X <- arima.sim(model = list(ar = c(1.2, -0.7), ma = 0.63), n = n,
#'  rand.gen = rnorm, n.start = n.start) + 7.7
#'
#'# Application of normCast()
#'result <- normCast(X = X, p = 2, q = 1, include.mean = TRUE, h = 5,
#'  plot = TRUE, xlim = c(1971, 2005), col = "deepskyblue4",
#'  type = "b", lty = 3, pch = 16, main = "Exemplary title")
#'result


normCast <- function(X, p = NULL, q = NULL, include.mean = FALSE, h = 1,
  alpha = 0.95, plot = FALSE, ...) {

  if (length(X) <= 1 || !all(!is.na(X)) || !is.numeric(X)) {
    stop("The argument 'X' must be a numeric vector with length ",
         "significantly > 1 and without NAs.")
  }
  if (length(h) != 1 || is.na(h) || !is.numeric(h) || h < 1) {
    stop("Argument 'h' must be a positive integer or at least a value >= 1.")
  }
  h <- floor(h)
  if (!(is.null(p) || (all(!is.na(p)) && is.numeric(p))) ||
      (is.numeric(p) && (length(p) > 1 || p < 0))) {
    stop("The argument 'p' must be an integer >= 0 or NULL.")
  }
  if (!(is.null(q) || (all(!is.na(q)) && is.numeric(q))) ||
      (is.numeric(q) && (length(q) > 1 || q < 0))) {
    stop("The argument 'q' must be an integer >= 0 or NULL.")
  }
  if (length(include.mean) != 1 || is.na(include.mean) ||
      !include.mean %in% c(TRUE, FALSE)) {
    stop("The argument 'include.mean' must be a single logical value ",
         "(TRUE or FALSE).")
  }
  if (length(alpha) != 1 || is.na(alpha) || !is.numeric(alpha) ||
      (alpha <= 0 || alpha >= 1)) {
    stop("The argument 'alpha' must be 0 < alpha < 1.")
  }
  if (length(plot) != 1 || is.na(plot) || !plot %in% c(TRUE, FALSE)) {
    stop("The argument 'plot' must be a single logical value (TRUE or FALSE).")
  }

  if (is.null(p) && is.null(q)) {
    message("Model selection in progress.")
    p.max <- 5
    q.max <- 5
    BICmat <- unname(critMatrix(X, p.max = p.max, q.max = q.max,
      include.mean = include.mean, criterion = "bic"))
    pq.opt <- which(BICmat == min(BICmat), arr.ind = TRUE) - 1
    p <- pq.opt[[1]]
    q <- pq.opt[[2]]
    message("Orders p=", p, " and q=", q, " were selected.")
  }
  if (is.null(p)) p <- 0
  if (is.null(q)) q <- 0
  p <- floor(p)
  q <- floor(q)

  suppressWarnings(pilot.est <- arima(X, order = c(p, 0, q),
    include.mean = include.mean))
  coefs <- pilot.est$coef
  lc <- length(coefs)
  if (include.mean == TRUE) {
    mu <- coefs[["intercept"]]
  } else {
    mu <- 0
  }
  if (p > 0) {
    ar <- head(coefs, p)
  } else {
    ar <- 0
  }
  if (q > 0) {
    ma <- head(coefs[(p + 1):lc], q)
  } else {
    ma <- 0
  }

  innov <- pilot.est$residuals
  sig2 <- pilot.est$sigma2

  X.fcst <- c(fcastCpp(X = X, innov = innov, ar = ar, ma = ma, mu = mu,
    h = h))
  c.coef <- maInfty(ar = ar, ma = ma, m = h - 1)
  sd.fcast <- sqrt(sig2 * cumsum(c.coef^2))
  alpha.s <- 1 - alpha
  q_norm <- qnorm(1 - alpha.s/2)
  m.out <- matrix(NA, ncol = h, nrow = 3)
  rname <- paste0(c(alpha.s/2, 1 - alpha.s/2) * 100, "%")
  rownames(m.out) <- c("fcast", rname)
  colnames(m.out) <- paste0("k=", 1:h)
  fi.low <- X.fcst - q_norm * sd.fcast
  fi.up <- X.fcst + q_norm * sd.fcast
  m.out[1, ] <- X.fcst
  m.out[2, ] <- fi.low
  m.out[3, ] <- fi.up

  if (plot == TRUE) {
    n <- length(X)
    dots <- list(...)
    if (is.null(dots[["x"]])) dots[["x"]] <- 1:n
    dots[["y"]] <- X
    tStep <- dots$x[2] - dots$x[1]
    x.fc <- seq(from = dots$x[n] + tStep, to = dots$x[n] + h * tStep,
      by = tStep)
    if (is.null(dots[["xlab"]])) dots[["xlab"]] <- "Time"
    if (is.null(dots[["ylab"]])) dots[["ylab"]] <- "Series"
    if (is.null(dots[["main"]])) dots[["main"]] <- "Forecasting results"
    if (is.null(dots[["type"]])) dots[["type"]] <- "l"
    if (is.null(dots[["xlim"]])) {
      xlim1 <- tail(dots$x, 6 * h)[1]
      dots[["xlim"]] <- c(xlim1, dots$x[n] + h * tStep)
    }
    if (is.null(dots[["ylim"]])) {
      x.all <- c(dots$x, x.fc)
      low.bound <- c(X, m.out[2, ])
      up.bound <- c(X, m.out[3, ])
      x.low <- sum(x.all < dots$xlim[1]) + 1
      x.up <- sum(x.all <= dots$xlim[2])
      dots[["ylim"]] <- c(min(low.bound[x.low:x.up]),
        max(up.bound[x.low:x.up]))
    }
    do.call(what = graphics::plot, args = dots)
    if (h >= 2) {
      polygon(c(x.fc, rev(x.fc)),
        c(m.out[2, ], rev(m.out[3, ])), col = col.alpha("gray", alpha = 0.8),
        border = NA)
      lines(x.fc, X.fcst, col = "red", lty = 1)
    } else {
      lines(c(x.fc, x.fc), c(m.out[2, ], m.out[3, ]),
        col = col.alpha("gray", alpha = 0.8), lty = 1)
      points(x.fc, m.out[1, ], col = "red", pch = 20)
    }
  }

  m.out
}
