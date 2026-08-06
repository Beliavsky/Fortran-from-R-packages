#' Multi-step-ahead regime and correlation forecasts
#'
#' Produces genuine \eqn{h}-step-ahead forecasts of the regime distribution and
#' the implied correlation matrix from the end of the estimation sample. Starting
#' from the terminal filtered regime probabilities \eqn{\xi_T}, the regime
#' distribution is propagated through the Markov chain,
#' \eqn{\pi_{T+k} = \pi_{T+k-1} P_{T+k}}, and the forecast correlations are the
#' regime-probability-weighted averages of the regime correlation matrices,
#' \eqn{\hat\rho_{T+k} = \sum_s \pi_{T+k}(s)\,\rho^{(s)}}.
#'
#' For \code{"noX"}/\code{"const"} the transition matrix is constant, so
#' \eqn{\pi_{T+k} = \xi_T P^k}. For \code{"tvtp"} each step uses a covariate row:
#' supply future covariates in \code{X_future} (an \eqn{h \times p} matrix); if it
#' is \code{NULL} the last in-sample covariate row is held constant and a message
#' is emitted.
#'
#' @param object An \code{"rsdc_fit"} object from \code{\link{rsdc_estimate}}.
#' @param horizon Integer forecast horizon \eqn{h \ge 1} (default \code{10}).
#' @param residuals Optional \eqn{T \times K} residuals used to obtain the
#'   terminal filtered state; defaults to the matrix stored on the fit.
#' @param X Optional in-sample covariate matrix (\code{"tvtp"}); defaults to the
#'   fit's stored matrix.
#' @param X_future Optional \eqn{h \times p} matrix of future covariates
#'   (\code{"tvtp"} only). If \code{NULL}, the last observed covariate row is
#'   reused at every horizon.
#'
#' @returns A list with:
#' \describe{
#'   \item{horizon}{Integer vector \code{1:h}.}
#'   \item{regime_probs}{\eqn{h \times N} matrix of forecast regime probabilities.}
#'   \item{predicted_correlations}{\eqn{h \times C} matrix (\eqn{C = K(K-1)/2}) of
#'     forecast correlations, one column per asset pair.}
#' }
#'
#' @seealso \code{\link{rsdc_forecast}} for the in-sample / 70-30 path.
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
#' rsdc_forecast_ahead(fit, horizon = 5)
#' }
#' @export
rsdc_forecast_ahead <- function(object, horizon = 10L, residuals = NULL,
                                X = NULL, X_future = NULL) {
  if (!inherits(object, "rsdc_fit")) stop("object must be an 'rsdc_fit'.")
  if (length(horizon) != 1L || !is.numeric(horizon) || !is.finite(horizon) ||
      horizon < 1 || horizon != as.integer(horizon))
    stop("horizon must be a single finite positive whole number.")
  horizon <- as.integer(horizon)
  method <- object$method; N <- object$N; K <- object$K
  C <- K * (K - 1) / 2
  if (is.null(residuals)) residuals <- object$residuals
  if (is.null(residuals)) stop("No residuals available; pass 'residuals'.")
  if (is.null(X) && method == "tvtp") X <- object$X
  rho_matrix <- object$correlations

  # Terminal filtered regime distribution xi_T.
  if (method == "const") {
    xi_T <- 1
  } else {
    hh <- rsdc_hamilton(residuals, if (method == "tvtp") X else NULL,
                        object$beta, rho_matrix, K, N,
                        if (method == "noX") object$transition_matrix else NULL)
    fp <- hh$filtered_probs
    xi_T <- if (!is.null(fp)) fp[, ncol(fp)] else rep(1 / N, N)
  }

  # Future transition matrices.
  P_at_k <- function(k) {
    if (method == "tvtp") {
      x <- if (!is.null(X_future)) X_future[k, ] else X[nrow(X), ]
      .rsdc_P_at(object$beta, x, N)
    } else if (method == "noX") {
      object$transition_matrix
    } else {
      matrix(1, 1, 1)
    }
  }
  if (method == "tvtp" && is.null(X_future))
    message("X_future not supplied; holding the last observed covariate row constant ",
            "across the forecast horizon.")
  if (method == "tvtp" && !is.null(X_future)) {
    X_future <- as.matrix(X_future)
    if (ncol(X_future) != ncol(X))
      stop(sprintf("X_future must have %d column(s) to match the covariate matrix X (got %d).",
                   ncol(X), ncol(X_future)))
    if (nrow(X_future) < horizon)
      stop("X_future must have at least 'horizon' rows.")
    # A non-finite covariate propagates through the softmax into NA transition
    # matrices and NA forecasts, with no error anywhere along the way.
    if (any(!is.finite(X_future[seq_len(horizon), , drop = FALSE])))
      stop("X_future must not contain NA/NaN/Inf values in the first ",
           horizon, " row(s).")
  }

  reg <- matrix(NA_real_, horizon, N)
  pc  <- matrix(NA_real_, horizon, C)
  pi_k <- as.numeric(xi_T)
  for (k in seq_len(horizon)) {
    pi_k <- as.numeric(pi_k %*% P_at_k(k))         # predictive regime distribution
    reg[k, ] <- pi_k
    pc[k, ]  <- as.numeric(crossprod(pi_k, rho_matrix))  # weighted regime correlations
  }
  list(horizon = seq_len(horizon),
       regime_probs = reg,
       predicted_correlations = pc)
}
