# broom-style tidiers and a ggplot2 autoplot method for rsdc_fit (added in 1.5-0).

# Non-standard evaluation column names used inside ggplot2::aes() below.
utils::globalVariables(c("time", "prob", "regime"))

#' Tidy a fitted RSDC model
#'
#' \pkg{broom}-style tidier returning a one-row-per-parameter coefficient table.
#'
#' @param x An \code{"rsdc_fit"} object.
#' @param ... Unused; for generic compatibility.
#' @returns A data frame with columns \code{term}, \code{estimate},
#'   \code{std.error}, \code{statistic} and \code{p.value}.
#' @examples
#' \donttest{
#' # Two persistent regimes: low (0.1) vs high (0.8) correlation
#' sim <- rsdc_simulate(n = 500, X = matrix(1, 500, 1),
#'                      beta = matrix(qlogis(0.9), 2, 1),
#'                      mu = matrix(0, 2, 2),
#'                      sigma = array(c(1, 0.1, 0.1, 1,
#'                                      1, 0.8, 0.8, 1), c(2, 2, 2)),
#'                      N = 2, seed = 2)
#' fit <- rsdc_estimate("noX", residuals = sim$observations, N = 2)
#' generics::tidy(fit)
#' }
#' @importFrom stats coef pnorm
#' @importFrom generics tidy glance augment
#' @exportS3Method generics::tidy
tidy.rsdc_fit <- function(x, ...) {
  est <- stats::coef(x)
  V   <- tryCatch(stats::vcov(x), error = function(e) NULL)
  se  <- if (!is.null(V)) .rsdc_safe_sqrt(diag(V)) else rep(NA_real_, length(est))
  stat <- est / se
  data.frame(term = names(est), estimate = as.numeric(est),
             std.error = as.numeric(se), statistic = as.numeric(stat),
             p.value = 2 * stats::pnorm(-abs(stat)), row.names = NULL)
}

#' Glance at a fitted RSDC model
#'
#' \pkg{broom}-style one-row model summary.
#'
#' @param x An \code{"rsdc_fit"} object.
#' @param ... Unused; for generic compatibility.
#' @returns A one-row data frame with \code{method}, \code{n_regimes},
#'   \code{logLik}, \code{AIC}, \code{BIC}, \code{df} and \code{nobs}.
#' @examples
#' \donttest{
#' # Two persistent regimes: low (0.1) vs high (0.8) correlation
#' sim <- rsdc_simulate(n = 500, X = matrix(1, 500, 1),
#'                      beta = matrix(qlogis(0.9), 2, 1),
#'                      mu = matrix(0, 2, 2),
#'                      sigma = array(c(1, 0.1, 0.1, 1,
#'                                      1, 0.8, 0.8, 1), c(2, 2, 2)),
#'                      N = 2, seed = 2)
#' fit <- rsdc_estimate("noX", residuals = sim$observations, N = 2)
#' generics::glance(fit)
#' }
#' @importFrom stats logLik AIC BIC nobs
#' @exportS3Method generics::glance
glance.rsdc_fit <- function(x, ...) {
  data.frame(method = x$method, n_regimes = x$N,
             logLik = as.numeric(stats::logLik(x)),
             AIC = stats::AIC(x), BIC = stats::BIC(x),
             df = x$npar, nobs = stats::nobs(x), row.names = NULL)
}

#' Augment data with fitted RSDC regime information
#'
#' \pkg{broom}-style augmentation: returns the estimation residuals together with
#' the smoothed regime probabilities and the most likely (Viterbi) regime path.
#'
#' @param x An \code{"rsdc_fit"} object (must carry its estimation residuals).
#' @param ... Unused; for generic compatibility.
#' @returns A data frame with one row per observation: the residual columns
#'   (\code{.resid1}, ...), smoothed regime probabilities
#'   (\code{.smoothed_p1}, ...) and the Viterbi state \code{.state}.
#' @examples
#' \donttest{
#' # Two persistent regimes: low (0.1) vs high (0.8) correlation
#' sim <- rsdc_simulate(n = 500, X = matrix(1, 500, 1),
#'                      beta = matrix(qlogis(0.9), 2, 1),
#'                      mu = matrix(0, 2, 2),
#'                      sigma = array(c(1, 0.1, 0.1, 1,
#'                                      1, 0.8, 0.8, 1), c(2, 2, 2)),
#'                      N = 2, seed = 2)
#' fit <- rsdc_estimate("noX", residuals = sim$observations, N = 2)
#' head(generics::augment(fit))
#' }
#' @exportS3Method generics::augment
augment.rsdc_fit <- function(x, ...) {
  res <- x$residuals
  if (is.null(res)) stop("augment() needs the estimation residuals on the fit.")
  df <- as.data.frame(res)
  names(df) <- paste0(".resid", seq_len(ncol(res)))
  if (x$method != "const" && x$N > 1L) {
    hh <- rsdc_hamilton(res, if (x$method == "tvtp") x$X else NULL,
                        x$beta, x$correlations, x$K, x$N,
                        if (x$method == "noX") x$transition_matrix else NULL)
    if (!is.null(hh$smoothed_probs)) {
      spt <- t(hh$smoothed_probs)
      colnames(spt) <- paste0(".smoothed_p", seq_len(ncol(spt)))
      df <- cbind(df, spt)
    }
    df$.state <- rsdc_viterbi(x, res, if (x$method == "tvtp") x$X else NULL)
  } else {
    df$.smoothed_p1 <- 1
    df$.state <- 1L
  }
  df
}

#' Plot a fitted RSDC model with \pkg{ggplot2}
#'
#' Plots the smoothed regime probabilities over the sample as stacked areas.
#' Requires the suggested \pkg{ggplot2} package.
#'
#' @param object An \code{"rsdc_fit"} object.
#' @param ... Unused; for generic compatibility.
#' @returns A \code{ggplot} object.
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   # Two persistent regimes: low (0.1) vs high (0.8) correlation
#'   sim <- rsdc_simulate(n = 500, X = matrix(1, 500, 1),
#'                        beta = matrix(qlogis(0.9), 2, 1),
#'                        mu = matrix(0, 2, 2),
#'                        sigma = array(c(1, 0.1, 0.1, 1,
#'                                        1, 0.8, 0.8, 1), c(2, 2, 2)),
#'                        N = 2, seed = 2)
#'   fit <- rsdc_estimate("noX", residuals = sim$observations, N = 2)
#'   ggplot2::autoplot(fit)
#' }
#' }
#' @exportS3Method ggplot2::autoplot
autoplot.rsdc_fit <- function(object, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for autoplot(); install it first.")
  res <- object$residuals
  if (is.null(res)) stop("autoplot() needs the estimation residuals on the fit.")
  if (object$method == "const" || object$N == 1L)
    stop("autoplot() is for switching models (N >= 2); 'const' has a single regime.")
  hh <- rsdc_hamilton(res, if (object$method == "tvtp") object$X else NULL,
                      object$beta, object$correlations, object$K, object$N,
                      if (object$method == "noX") object$transition_matrix else NULL)
  sp <- hh$smoothed_probs
  if (is.null(sp)) stop("Hamilton filter returned no smoothed probabilities.")
  N <- nrow(sp); Tn <- ncol(sp)
  long <- data.frame(
    time   = rep(seq_len(Tn), times = N),
    regime = factor(rep(seq_len(N), each = Tn)),
    prob   = as.numeric(t(sp)))
  ggplot2::ggplot(long, ggplot2::aes(x = time, y = prob, fill = regime)) +
    ggplot2::geom_area(position = "stack") +
    ggplot2::labs(x = "Time", y = "Smoothed probability", fill = "Regime") +
    ggplot2::theme_minimal()
}
