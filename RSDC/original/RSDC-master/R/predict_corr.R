#' Forecast Covariance/Correlation Paths from an RSDC Model
#'
#' Generates per-period correlation and covariance matrices from a fitted model:
#' \code{"const"} (constant correlation), \code{"noX"} (fixed transition matrix), or
#' \code{"tvtp"} (time-varying transition probabilities).
#'
#' @param method Character. One of \code{"tvtp"}, \code{"noX"}, \code{"const"}.
#' @param N Integer. Number of regimes (ignored for \code{"const"}).
#' @param residuals Numeric matrix \eqn{T \times K} used to compute correlations or run the filter.
#' @param X Optional numeric matrix \eqn{T \times p} (required for \code{"tvtp"}).
#' @param final_params List with fitted parameters (e.g., from \code{\link{rsdc_estimate}}):
#'   must include \code{correlations}, and either \code{transition_matrix} (\code{"noX"})
#'   or \code{beta} (\code{"tvtp"}); include \code{log_likelihood} for BIC computation.
#' @param sigma_matrix Numeric matrix \eqn{T \times K} of forecasted standard deviations.
#' @param value_cols Character/integer vector of columns in \code{sigma_matrix} that define asset order.
#' @param out_of_sample Logical. If \code{TRUE}, use a fixed 70/30 split; otherwise use the whole sample.
#' @param control Optional list; supports \code{threshold} (in (0,1), default \code{0.7}).
#'
#' @return
#' \describe{
#'   \item{\code{smoothed_probs}}{\eqn{N \times T^\ast} smoothed probabilities (\code{"noX"}/\code{"tvtp"} only).}
#'   \item{\code{sigma_matrix}}{\eqn{T^\ast \times K} slice aligned to the forecast horizon.}
#'   \item{\code{cov_matrices}}{List of \eqn{K \times K} covariance matrices \eqn{\Sigma_t = D_t R_t D_t}.}
#'   \item{\code{predicted_correlations}}{\eqn{T^\ast \times \binom{K}{2}} pairwise correlations in \code{combn(K, 2)} order.}
#'   \item{\code{BIC}}{Information criterion \eqn{\log(n)\,k - 2\,\ell}.
#'     When \code{out_of_sample = TRUE} and \code{final_params$log_likelihood_oos} is non-\code{NULL},
#'     uses the OOS log-likelihood and OOS sample size; otherwise uses the IS values.
#'     \strong{Caveat:} with \code{out_of_sample = TRUE} this is an out-of-sample predictive
#'     score (held-out log-likelihood penalized by the IS parameter count), \emph{not} a
#'     conventional in-sample BIC; use it only to compare models fit on the same split, and
#'     do not interpret it as a textbook BIC. The backends (\code{\link{rsdc_estimate}})
#'     hardcode a 70/30 split, so the OOS log-likelihood is only fully consistent with
#'     \code{rsdc_forecast} when the default \code{threshold = 0.7} is used.}
#'   \item{\code{y}}{\eqn{T^\ast \times K} residual matrix aligned to the forecast horizon.}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{Forecast horizon:} If \code{out_of_sample = TRUE}, filter on the first \code{threshold}
#'         fraction and forecast on the remainder.
#'   \item \strong{Correlation paths:}
#'         \itemize{
#'           \item \code{"const"} — constant correlation from \code{final_params$correlations[1,]}, repeated across time.
#'           \item \code{"noX"}/\code{"tvtp"} — regime-probability weighted average of regime
#'                 correlations. Out-of-sample (\code{out_of_sample = TRUE}) the weights are the
#'                 \emph{filtered} probabilities \eqn{\Pr(S_t \mid \Omega_t)}, which use no
#'                 \emph{future} information beyond time \eqn{t} (the OOS window is initialized from
#'                 the in-sample terminal filtered state). Note these condition on the
#'                 \emph{contemporaneous} \eqn{U_t}, so they are a real-time \emph{nowcast}, not a
#'                 one-step-ahead forecast formed before the covariate is seen; for that use
#'                 \code{\link{rsdc_forecast_ahead}}. The full-sample case uses \emph{smoothed}
#'                 probabilities \eqn{\Pr(S_t \mid \Omega_T)}. The returned
#'                 \code{smoothed_probs} field always contains the smoothed probabilities.
#'         }
#'   \item \strong{Covariance build:} Reconstruct \eqn{R_t} from the pairwise vector (columns ordered by
#'         \code{combn(K, 2)}), set \eqn{D_t = \mathrm{diag}(\sigma_{t,1},\dots,\sigma_{t,K})}, and
#'         \eqn{\Sigma_t = D_t R_t D_t}.
#'   \item \strong{BIC:} Parameter count \eqn{k} is
#'         \code{N * ncol(X) + N * K * (K - 1) / 2} for \code{"tvtp"} with \eqn{N=2},
#'         \code{N * (N-1) * ncol(X) + N * K * (K - 1) / 2} for \code{"tvtp"} with \eqn{N \ge 3},
#'         \code{N * (N - 1) + N * K * (K - 1) / 2} for \code{"noX"},
#'         and \code{K * (K - 1) / 2} for \code{"const"}.
#' }
#'
#' @examples
#' set.seed(123)
#' T <- 60; K <- 3; N <- 2
#' y <- scale(matrix(rnorm(T*K), T, K))
#' vols <- matrix(0.2 + 0.05*abs(sin(seq_len(T)/7)), T, K)
#' rho <- rbind(c(0.10, 0.05, 0.00), c(0.60, 0.40, 0.30))
#' Pfix <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
#' rsdc_forecast("noX", N, y, NULL,
#'               list(correlations = rho, transition_matrix = Pfix, log_likelihood = -200),
#'               vols, 1:K)
#' @seealso \code{\link{rsdc_hamilton}}, \code{\link{rsdc_estimate}},
#'   \code{\link{rsdc_minvar}}, \code{\link{rsdc_maxdiv}}
#'
#' @importFrom stats cor
#' @importFrom utils combn
#' @export
rsdc_forecast <- function(method = c("tvtp", "noX", "const"),
                         N,
                         residuals,
                         X = NULL,
                         final_params,
                         sigma_matrix,
                         value_cols,
                         out_of_sample = FALSE,
                         control = list()) {

  con <- list(threshold = 0.7)
  con[names(control)] <- control
  thresh <- con$threshold
  if (length(thresh) != 1L || !is.numeric(thresh) || !is.finite(thresh) ||
      thresh <= 0 || thresh >= 1)
    stop("control$threshold must be a single number strictly between 0 and 1.")

  method <- match.arg(method)
  if (method == "tvtp" && is.null(X))
    stop("X must be provided for method = 'tvtp'.")
  K <- ncol(residuals)
  n_obs <- nrow(residuals)

  # Input dimension validation (mirrors rsdc_minvar/rsdc_maxdiv)
  if (!is.matrix(residuals)) stop("residuals must be a numeric matrix.")
  if (length(value_cols) != K)
    stop(sprintf("length(value_cols) = %d must equal K = ncol(residuals) = %d.",
                 length(value_cols), K))
  if (nrow(sigma_matrix) != n_obs)
    stop(sprintf("nrow(sigma_matrix) = %d must equal nrow(residuals) = %d.",
                 nrow(sigma_matrix), n_obs))
  if (ncol(sigma_matrix) < K)
    stop(sprintf("sigma_matrix must have at least K = %d columns (got %d).",
                 K, ncol(sigma_matrix)))
  if (!is.null(X) && nrow(X) != n_obs)
    stop(sprintf("nrow(X) = %d must equal nrow(residuals) = %d.", nrow(X), n_obs))

  if (out_of_sample) {
    is_cut <- round(thresh * n_obs)
    if (is_cut < 2L || n_obs - is_cut < 1L)
      stop("Sample too small for the in-sample/out-of-sample split.")
    residuals_is <- residuals[1:is_cut, , drop = FALSE]
    residuals_oos <- residuals[(is_cut + 1):n_obs, , drop = FALSE]
    sigma_matrix_oos <- sigma_matrix[(is_cut + 1):n_obs, , drop = FALSE]
    X_is <- if (!is.null(X)) X[1:is_cut, , drop = FALSE] else NULL
    X_oos <- if (!is.null(X)) X[(is_cut + 1):n_obs, , drop = FALSE] else NULL
  } else {
    residuals_oos <- residuals
    sigma_matrix_oos <- sigma_matrix
    X_oos <- X
  }

  # Default outputs
  smoothed_probs <- NULL
  predicted_correl <- NULL

  if (method == "const") {
    # Use the fitted constant correlation from the estimated model
    rho_vec <- as.numeric(final_params$correlations[1, ])
    predicted_correl <- matrix(rep(rho_vec, nrow(residuals_oos)),
                               ncol = K * (K - 1) / 2, byrow = TRUE)
  } else {
    # Use Hamilton filter to get smoothed_probs
    if (out_of_sample) {
      # Pass 1: filter IS data to get the terminal filtered state
      result_is <- rsdc_hamilton(
        y          = residuals_is,
        X          = if (method == "tvtp") X_is else NULL,
        beta       = if (method == "tvtp") final_params$beta else NULL,
        rho_matrix = final_params$correlations,
        K = K, N = N,
        P          = if (method == "noX") final_params$transition_matrix else NULL
      )
      xi_terminal <- if (!is.null(result_is$filtered_probs))
        result_is$filtered_probs[, ncol(result_is$filtered_probs)]
      else rep(1 / N, N)

      # Pass 2: filter OOS data initialized from the IS terminal state
      result_filter <- rsdc_hamilton(
        y          = residuals_oos,
        X          = if (method == "tvtp") X_oos else NULL,
        beta       = if (method == "tvtp") final_params$beta else NULL,
        rho_matrix = final_params$correlations,
        K = K, N = N,
        P          = if (method == "noX") final_params$transition_matrix else NULL,
        xi_init    = xi_terminal
      )
    } else {
      # Full-sample: filter on the complete dataset
      result_filter <- rsdc_hamilton(
        y          = residuals,
        X          = if (method == "tvtp") X else NULL,
        beta       = if (method == "tvtp") final_params$beta else NULL,
        rho_matrix = final_params$correlations,
        K = K, N = N,
        P          = if (method == "noX") final_params$transition_matrix else NULL
      )
    }

    # Genuine OOS forecasts must not use future data: smoothed probs condition on the
    # entire OOS window (look-ahead). Use FILTERED probs out-of-sample; keep smoothed
    # for in-sample (full-sample) reconstruction.
    weight_probs   <- if (out_of_sample) result_filter$filtered_probs
                      else               result_filter$smoothed_probs
    smoothed_probs <- result_filter$smoothed_probs
    if (is.null(weight_probs))
      stop("Hamilton filter returned -Inf likelihood; check that rho_matrix is positive definite.")

    # Weighted correlation forecast — size derived from filter output, not from residuals_oos
    cor_pairs <- combn(K, 2, simplify = FALSE)
    predicted_correl <- matrix(0, nrow = ncol(weight_probs), ncol = length(cor_pairs))

    for (t in 1:nrow(predicted_correl)) {
      for (s in 1:N) {
        predicted_correl[t, ] <- predicted_correl[t, ] +
          weight_probs[s, t] * final_params$correlations[s, ]
      }
    }
  }

  # Covariance matrix construction
  cor_pairs <- combn(K, 2, simplify = FALSE)
  cov_matrices <- lapply(1:nrow(sigma_matrix_oos), function(t) {
    sigma_vals <- as.numeric(sigma_matrix_oos[t, value_cols])
    R <- diag(K)

    for (i in seq_along(cor_pairs)) {
      pair <- cor_pairs[[i]]
      R[pair[1], pair[2]] <- R[pair[2], pair[1]] <- predicted_correl[t, i]
    }

    # Same validation as the portfolio optimizers: a negative or non-finite
    # "standard deviation" flips covariance signs through D %*% R %*% D, and an
    # indefinite R is not a correlation matrix. Fail at the offending period
    # rather than returning a silently wrong covariance.
    .rsdc_check_cov_inputs(sigma_vals, R, t)
    D <- diag(sigma_vals)
    D %*% R %*% D
  })

  # BIC calculation
  k <- switch(
    method,
    "tvtp" = (if (N == 2) N * ncol(X) else N * (N - 1L) * ncol(X)) +
             N * K * (K - 1) / 2,
    "noX"  = N * (N - 1) + N * K * (K - 1) / 2,
    "const" = K * (K - 1) / 2
  )
  if (out_of_sample && !is.null(final_params$log_likelihood_oos)) {
    n_bic  <- nrow(residuals_oos)
    ll_bic <- final_params$log_likelihood_oos
  } else {
    n_bic  <- if (out_of_sample) is_cut else n_obs
    ll_bic <- final_params$log_likelihood
  }
  if (is.null(ll_bic)) {
    # `2 * NULL` would yield numeric(0); return a scalar NA so downstream code
    # expecting a length-1 BIC does not break silently.
    warning("BIC is NA: log_likelihood missing from final_params.")
    ll_bic <- NA_real_
  }
  BIC <- log(n_bic) * k - 2 * ll_bic

  list(
    smoothed_probs = smoothed_probs,
    sigma_matrix = sigma_matrix_oos,
    cov_matrices = cov_matrices,
    predicted_correlations = predicted_correl,
    BIC = BIC,
    y = residuals_oos
  )
}
