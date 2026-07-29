#' Forecasting Function for ARMA Models via Bootstrap
#'
#' Point forecasts and the respective forecasting intervals for
#' autoregressive-moving-average (ARMA) models can be calculated, the latter
#' via bootstrap, by means of this function.
#'
#'@param X a numeric vector that contains the time series that is assumed to
#'follow an ARMA model ordered from past to present.
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
#'@param include.mean a logical value; if set to \code{TRUE}, the mean of the
#'series is also also estimated; if set to \code{FALSE}, \eqn{E(X_t) = 0} is
#'assumed; is set to \code{FALSE} by default.
#'@param n.start an integer that defines the 'burn-in' number
#'of observations for the simulated ARMA series via bootstrap; is set to
#'\code{1000} by default; decimal numbers will be rounded off to integers.
#'@param h an integer that represents the forecasting horizon; if \eqn{n} is
#'the number of observations, point forecasts and forecasting intervals will be
#'obtained for the time points \eqn{n + 1} to \eqn{n + h}; is set to
#'\code{h = 1} by default; decimal numbers will be rounded off to integers.
#'@param it an integer that represents the total number of iterations, i.e.,
#'the number of simulated series; is set to \code{10000} by default; decimal
#'numbers will be rounded off to integers.
#'@param pb a logical value; for \code{pb = TRUE}, a progress bar will be shown
#'in the console.
#'@param cores an integer value >0 that states the number of (logical) cores to
#'use in the bootstrap (or \code{NULL}); the default is the maximum number of
#'available cores
#'(via \code{\link[future:availableCores]{future::availableCores}}); for
#'\code{cores = NULL}, parallel computation is disabled.
#'@param alpha a numeric vector of length 1 with \eqn{0 < } \code{alpha}
#'\eqn{ < 1}; the forecasting intervals will be obtained based on the
#'confidence level (\eqn{100}\code{alpha})-percent; is set to
#'\code{alpha = 0.95} by default, i.e., a \eqn{95}-percent confidence level.
#'@param export.error a single logical value; if the argument is set to
#'\code{TRUE}, a list is returned instead of a matrix (\code{FALSE}); the
#'first element of the list is the usual forecasting matrix, whereas the second
#'element is a matrix with \code{h} columns, where each column represents
#'the calculated forecasting errors for the respective future time point
#'\eqn{n + 1, n + 2, ..., n + h}; is set to \code{FALSE} by default.
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
#'version 1.1.0. For a given time series \eqn{X_t}{X_[t]}, \eqn{t = 1, 2, ..., n},
#'the point forecasts and the respective forecasting intervals will be
#'calculated. It is assumed that the series follows an ARMA(\eqn{p,q}) model
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
#'The forecasting intervals on the other hand are obtained by a forward
#'bootstrap method that was introduced by Pan and Politis (2016) for
#'autoregressive models and extended by Lu and Wang (2020) for applications to
#'autoregressive-moving-average models.
#'For this purpose, let \eqn{l} be the number of the current bootstrap
#'iteration. Based on the demeaned residuals of the initial ARMA estimation,
#'different innovation series \eqn{\epsilon_{l,t}^{s}}{\epsilon^[s]_[l,t]} will
#'be sampled.  The initial coefficient estimates and the sampled innovation
#'series are then used to simulate a variety of series
#'\eqn{X_{l,t}^{s}}{X^[s]_[l,t]}, from which again coefficient estimates will
#'be obtained. With these newly obtained estimates, proxy residual series
#'\eqn{\hat{\epsilon}_{l,t}^{s}}{hat[\epsilon]^[s]_[l,t]} are calculated for
#'the original series \eqn{X_t}{X_[t]}. Subsequently, point forecasts for the
#'time points \eqn{n + 1} to \eqn{n + h} are obtained for each iteration
#'\eqn{l} based on the original series \eqn{X_t}{X_[t]}, the newly obtained
#'coefficient forecasts and the proxy residual series
#'\eqn{\epsilon_{l,t}^{s}}{hat[\epsilon]^[s]_[l,t]}.
#'Simultaneously, "true" forecasts, i.e., true future observations, are
#'simulated. Within each iteration, the difference between the simulated true
#'forecast and the bootstrapped point forecast is calculated and saved for each
#'future time point \eqn{n + 1} to \eqn{n + h}. The result for these time
#'points are simulated empirical values of the forecasting error. Denote by
#'\eqn{q_k(.)}{q_[k](.)} the quantile of the empirical distribution for the
#'future time point \eqn{n + k}. Given a predefined confidence level
#'\code{alpha}, define \eqn{\alpha_s = (1 -} \code{alpha}\eqn{)/2}. The
#'bootstrapped forecasting interval is then
#'\deqn{[\hat{X}_{n + k} + q_k(\alpha_s), \hat{X}_{n + k} + q_k(1 -
#'\alpha_s)],}{[hat[X]_[n+k] + q_[k](alpha_[s]), hat[X]_[n+k] + q_[k](1 -
#'alpha_[s])],}
#'i.e., the forecasting intervals are given by the sum of the respective point
#'forecasts and quantiles of the respective bootstrapped forecasting error
#'distributions.
#'
#'The function \code{bootCast} allows for different adjustments to
#'the forecasting progress. At first, a vector with the values of the observed
#'time series ordered from past to present has to be passed to the argument
#'\code{X}. Orders \eqn{p} and \eqn{q} of the underlying ARMA process can be
#'defined via the arguments \code{p} and \code{q}. If only one of these orders
#'is inserted by the user, the other order is automatically set to \code{0}. If
#'none of these arguments are defined, the function will choose orders based on
#'the Bayesian Information Criterion (BIC) for
#'\eqn{0 \leq p,q \leq 5}{0 \le p,q \le 5}. Via the logical argument
#'\code{include.mean} the user can decide, whether to consider the mean of the
#'series within the estimation process. By means of \code{n.start}, the number
#'of "burn-in" observations for the simulated ARMA processes can be regulated.
#'These observations are usually used for the processes to build up and then
#'omitted. Furthermore, the argument \code{h} allows for the definition of the
#'maximum future time point \eqn{n + h}. Point forecasts and forecasting
#'intervals will be returned for the time points \eqn{n + 1} to \eqn{n + h}.
#'\code{it} corresponds to the number of bootstrap iterations. We recommend a
#'sufficiently high number of repetitions for maximum accuracy of the results.
#'Another argument is \code{alpha}, which is the equivalent of the confidence
#'level considered within the calculation of the forecasting intervals, i.e.,
#'the quantiles \eqn{(1 - } \code{alpha}\eqn{)/2} and \eqn{1 - (1 - }
#'\code{alpha}\eqn{)/2} of the bootstrapped forecasting error distribution
#'will be obtained.
#'
#'Since this bootstrap approach needs a lot of computation time, especially for
#'series with high numbers of observations and when fitting models with many
#'parameters, parallel computation of the bootstrap iterations is enabled.
#'With \code{cores}, the number of cores can be defined with an integer.
#'Nonetheless, for \code{cores = NULL}, no cluster is created and therefore
#'the parallel computation is disabled. Note that the bootstrapped results are
#'fully reproducible for all cluster sizes. The progress of the bootstrap can
#'be observed in the R console, where a progress bar and the estimated
#'remaining time are displayed for \code{pb = TRUE}.
#'
#'If the argument \code{export.error} is set to \code{TRUE}, the output of
#'the function is a list instead of a matrix with additional information on
#'the simulated forecasting errors. For more information see the section
#'\emph{Value}.
#'
#'For simplicity, the function also incorporates the possibility to directly
#'create a plot of the output, if the argument \code{plot} is set to
#'\code{TRUE}. By the additional and optional arguments \code{...}, further
#'arguments of the standard plot function can be implemented to shape the
#'returned plot.
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
#'@return The function returns a \eqn{3} by \eqn{h} matrix with its columns
#'representing the future time points and the point forecasts, the lower
#'bounds of the forecasting intervals and the upper bounds of the
#'forecasting intervals as the rows. If the argument \code{plot} is set to
#'\code{TRUE}, a plot of the forecasting results is created.
#'
#'If \code{export.error = TRUE} is selected, a list with the following
#'elements is returned instead.
#'
#'\describe{
#'\item{fcast}{the \eqn{3} by \eqn{h} matrix forecasting matrix with point
#'forecasts and bounds of the forecasting intervals.}
#'\item{error}{a \code{it} by \eqn{h} matrix, where each column represents a
#'future time point \eqn{n + 1, n + 2, ..., n + h}; in each column the
#'respective \code{it} simulated forecasting errors are saved.}
#'}
#'
#'@export
#'
#'@importFrom graphics lines points polygon
#'@importFrom stats arima quantile arima.sim
#'@importFrom utils head tail
#'@importFrom future plan sequential multisession
#'@importFrom future.apply future_lapply
#'@importFrom progressr progressor with_progress handler_progress
#'@importFrom progress progress_bar
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
#' Lu, X., and Wang, L. (2020). Bootstrap prediction interval for ARMA models
#' with unknown orders. REVSTAT–Statistical Journal, 18:3, 375-396.
#'
#' Pan, L. and Politis, D. N. (2016). Bootstrap prediction intervals for linear,
#' nonlinear and nonparametric autoregressions. In: Journal of Statistical
#' Planning and Inference 177, pp. 1-27.
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
#'# Function for drawing from a demeaned chi-squared distribution
#'rchisq0 <- function(n, df, npc = 0) {
#'  rchisq(n, df, npc) - df
#'}
#'
#'# Simulation of the underlying process
#'n <- 2000
#'n.start = 1000
#'set.seed(23)
#'X <- arima.sim(model = list(ar = c(1.2, -0.7), ma = 0.63), n = n,
#'  rand.gen = rchisq0, n.start = n.start, df = 3) + 13.1
#'
#'# Quick application with low number of iterations
#'# (not recommended in practice)
#'result <- bootCast(X = X, p = 2, q = 1, include.mean = TRUE,
#'  n.start = n.start, h = 5, it = 10, cores = 2, plot = TRUE,
#'  lty = 3, col = "forestgreen", xlim = c(1950, 2005), type = "b",
#'  main = "Exemplary title", pch = "*")
#'result
#'
#'### Example 2: Application with more iterations ###
#'\dontrun{
#'result2 <- bootCast(X = X, p = 2, q = 1, include.mean = TRUE,
#'  n.start = n.start, h = 5, it = 10000, cores = 2, plot = TRUE,
#'  lty = 3, col = "forestgreen", xlim = c(1950, 2005),
#'  main = "Exemplary title")
#'result2
#'
#'}

bootCast <- function(X, p = NULL, q = NULL, include.mean = FALSE,
  n.start = 1000, h = 1, it = 10000, pb = TRUE,
  cores = future::availableCores(), alpha = 0.95,
  export.error = FALSE, plot = FALSE, ...) {

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
  if (length(n.start) != 1 || is.na(n.start) || !is.numeric(n.start) ||
      n.start < 0) {
    stop("The argument 'n.start' must be a single non-negative integer value.")
  }
  n.start <- floor(n.start)
  if (length(it) != 1 || is.na(it) || !is.numeric(it) || it < 1) {
    stop("The argument 'it' must be a single positive integer value.")
  }
  it <- floor(it)

  if (length(pb) != 1 || is.na(pb) || !is.logical(pb)) {
    stop("The argument 'pb' must be a single logical value.")
  }
  if (!(length(cores) %in% c(0, 1)) || !(is.null(cores) || !is.na(cores) || is.numeric(cores))) {
    stop("The argument 'cores' must be a single integer value.")
  }
  if (length(plot) != 1 || is.na(plot) || !plot %in% c(TRUE, FALSE)) {
    stop("The argument 'plot' must be a single logical value (TRUE or FALSE).")
  }
  if (length(export.error) != 1 || is.na(export.error) ||
    !export.error %in% c(TRUE, FALSE)) {
    stop("The argument 'export.error' must be a single logical value ",
         "(TRUE or FALSE).")
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

  bootSet <- bootSetup(X, p, q, include.mean, h, n.start)
  bootF <- bootSet[[1]]
  X.fcast <- bootSet[[2]]
  old.plan <- future::plan()
  if (is.null(cores)) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession, workers = cores)
  }
  on.exit(future::plan(old.plan))

  progressr::with_progress({
    if (pb) {
      prog <- progressr::progressor(it / 100)
    } else {
      prog <- function() {invisible(NULL)}
    }
    errs <- future.apply::future_lapply(1:it,
      function(x, ...) {
        out <- bootF(x)
        if (x %% 100 == 0) prog()
        out
      }, future.seed = TRUE)

  }, handlers = progressr::handler_progress)

  err.mat <- matrix(unlist(errs), nrow = it, ncol = h, byrow = TRUE)

  alpha.s <- 1 - alpha
  quants <- apply(err.mat, MARGIN = 2, quantile,
                  probs = c(alpha.s / 2, 1 - alpha.s / 2))
  for (i in 1:2) {
    quants[i, ] <- quants[i, ] + X.fcast
  }
  out <- rbind(fcast = X.fcast, quants)
  colnames(out) <- paste0("k=", 1:h)

  if (plot == TRUE) {
    n <- length(X)
    dots <- list(...)
    if (is.null(dots[["x"]])) dots[["x"]] <- 1:n
    dots[["y"]] <- X
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
      low.bound <- c(X, out[2, ])
      up.bound <- c(X, out[3, ])
      x.low <- sum(x.all < dots$xlim[1]) + 1
      x.up <- sum(x.all <= dots$xlim[2])
      dots[["ylim"]] <- c(min(low.bound[x.low:x.up]),
        max(up.bound[x.low:x.up]))
    }
    do.call(what = graphics::plot, args = dots)
    if (h >= 2) {
      polygon(c(x.fc, rev(x.fc)),
        c(out[2, ], rev(out[3, ])), col = col.alpha("gray", alpha = 0.8),
        border = NA)
      lines(x.fc, X.fcast, col = "red", lty = 1)
    } else {
      lines(c(x.fc, x.fc), c(out[2, ], out[3, ]),
        col = col.alpha("gray", alpha = 0.8))
      points(x.fc, out[1, ], col = "red", pch = 20)
    }
  }

  if (export.error == TRUE) {
    out <- list(fcast = out, error = err.mat)
  }

  out

}

bootSetup <- function(X, p, q, include.mean, h, n.start) {
  force(n.start)
  n <- length(X)
  estMain <- suppressWarnings(stats::arima(X, order = c(p, 0, q),
    include.mean = include.mean))
  innov <- estMain$residuals
  Fi <- innov - mean(innov)
  coefs <- estMain$coef
  rm(estMain)
  ar <- numeric(0)
  ma <- numeric(0)
  ar.t <- 0
  ma.t <- 0
  mu <- 0
  if (p > 0) {
    ar <- ar.t <- coefs[1:p]
  }
  if (q > 0) {
    ma <- ma.t <- coefs[(p + 1):(p + q)]
  }
  if (include.mean) {
    mu <- coefs[["intercept"]]
  }
  rm(coefs)
  X.fcast <- c(fcastCpp(X, innov, ar.t, ma.t, mu, h))

  ar.star <- 0
  ma.star <- 0
  mu.star <- 0

  list(
  function(x) {
    eps_boot <- sample(Fi, size = n.start + n + h, replace = TRUE)
    X.star <- as.numeric(stats::arima.sim(model = list(ar = ar, ma = ma),
      n = n, n.start = n.start, innov = eps_boot[(n.start + 1):(n.start + n)],
      start.innov = eps_boot[1:n.start])) + mu
    est <- suppressWarnings(
      tryCatch(
        stats::arima(X.star, order = c(p, 0, q), include.mean = include.mean),
        error = function(c1) {
          stats::arima(X.star, order = c(p, 0, q),
            include.mean = include.mean, method = "ML")
        }
      )
    )

    coef <- est$coef
    if (p > 0) {
      ar.star <- coef[1:p]
    }
    if (q > 0) {
      ma.star <- coef[(p + 1):(p + q)]
    }
    if (include.mean) {
      mu.star <- coef[["intercept"]]
    }

    eps.star <- suppressWarnings(stats::arima(X, order = c(p, 0, q),
      fixed = c(head(ar.star, p), head(ma.star, q),
        mu.star))[["residuals"]])
    X.star.hat <- c(fcastCpp(X, eps.star, ar.star, ma.star, mu.star, h))
    X.true <- c(tfcastCpp(X, innov,
      eps_boot[(n.start + n + 1):(n.start + n + h)], ar.t, ma.t, mu, h))
    X.true - X.star.hat
  },
  X.fcast)
}
