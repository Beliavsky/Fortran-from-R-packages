#' Uncertainty bands for the predicted correlation path
#'
#' Propagates parameter uncertainty into the (smoothed-probability-weighted)
#' predicted correlation path of a fitted \code{"rsdc_fit"} model. Parameter
#' vectors are drawn from the asymptotic sampling distribution
#' \eqn{\mathcal{N}(\hat\theta, \widehat{\mathrm{Var}}(\hat\theta))} (using the
#' covariance returned by \code{\link{vcov.rsdc_fit}}); for each draw the Hamilton
#' filter/smoother is rerun and the predicted correlation path recomputed. The
#' pointwise quantiles of these paths form the band.
#'
#' This reflects estimation (parameter) uncertainty in the correlation
#' parameters and transition dynamics; it does not add the irreducible
#' regime-classification uncertainty already summarized by the smoothed
#' probabilities.
#'
#' Draws are taken from the \emph{unconstrained} Gaussian approximation; draws
#' that yield an invalid model (e.g. \eqn{|\rho| \ge 1}, a non-positive-definite
#' correlation matrix, or an invalid transition matrix) are dropped with a
#' warning and the band is computed from the remaining draws. When the estimate
#' sits close to a parameter bound this truncates the sampling distribution, so
#' the reported band can be somewhat too narrow there; the attribute
#' \code{B_used} reports how many draws survived.
#'
#' @param object An \code{"rsdc_fit"} object (from \code{\link{rsdc_estimate}})
#'   estimated with \code{control = list(compute_se = TRUE)} so a covariance is
#'   available.
#' @param residuals Optional \eqn{T \times K} residual matrix; defaults to the
#'   matrix stored on the fit.
#' @param X Optional \eqn{T \times p} covariate matrix (\code{"tvtp"} only);
#'   defaults to the matrix stored on the fit.
#' @param B Integer. Number of parameter draws (default \code{500}).
#' @param level Confidence level for the pointwise band (default \code{0.95}).
#' @param seed Optional integer seed for reproducibility.
#' @param cores Integer (default 1). Number of cores used to rerun the filter
#'   over the draws in parallel (via the \pkg{parallel} package). All draws
#'   are generated up front, so the result is identical for any number of
#'   cores.
#'
#' @returns A list with one element per correlation pair
#'   (\eqn{C = K(K-1)/2} of them), each a \eqn{T \times 3} matrix with columns
#'   \code{fit}, \code{lower}, \code{upper}; plus attributes \code{level} and
#'   \code{B_used} (draws that yielded a valid path).
#'
#' @seealso \code{\link{rsdc_forecast}}, \code{\link{rsdc_bootstrap}}
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
#' bands <- rsdc_corr_bands(fit, B = 100, seed = 1)
#' head(bands[[1]])
#' }
#' @importFrom stats quantile
#' @export
rsdc_corr_bands <- function(object, residuals = NULL, X = NULL, B = 500L,
                            level = 0.95, seed = NULL, cores = 1) {
  if (!inherits(object, "rsdc_fit")) stop("object must be an 'rsdc_fit'.")
  if (!is.numeric(B) || length(B) != 1L || B < 2) stop("B must be a single integer >= 2.")
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("level must be a single number in (0, 1).")
  if (is.null(object$vcov))
    stop("rsdc_corr_bands needs a covariance matrix; refit with ",
         "control = list(compute_se = TRUE).")
  method <- object$method; N <- object$N; K <- object$K
  C <- K * (K - 1) / 2
  if (is.null(residuals)) residuals <- object$residuals
  if (is.null(residuals)) stop("No residuals available; pass 'residuals'.")
  if (is.null(X) && method == "tvtp") X <- object$X
  par <- object$par; V <- object$vcov
  if (!is.null(seed)) set.seed(seed)

  # Predicted correlation path (C x T) at a parameter vector, or NULL if invalid.
  path_of <- function(p) {
    up <- .rsdc_unpack_par(p, method, X, K, N)
    if (any(abs(up$rho_matrix) >= 1)) return(NULL)
    hh <- tryCatch(
      rsdc_hamilton(residuals, if (method == "tvtp") X else NULL,
                    up$beta, up$rho_matrix, K, N, up$P),
      error = function(e) NULL)
    sp <- hh$smoothed_probs
    if (is.null(sp)) return(NULL)
    pc <- matrix(0, C, ncol(sp))
    for (s in seq_len(N)) pc <- pc + outer(up$rho_matrix[s, ], sp[s, ])
    pc
  }

  point <- path_of(par)
  if (is.null(point)) stop("Could not compute the point-estimate correlation path.")
  Tn <- ncol(point)

  draws <- mvtnorm::rmvnorm(B, mean = par, sigma = V)
  paths <- .rsdc_lapply(seq_len(B), function(b) path_of(draws[b, ]),
                        cores = cores)
  paths <- paths[!vapply(paths, is.null, logical(1))]
  B_used <- length(paths)
  if (B_used < 2L)
    stop("Too few valid parameter draws produced a finite path; increase B or check the fit.")
  if (B_used < B)
    warning(sprintf(paste("%d of %d parameter draws produced an invalid (e.g. non-positive-definite)",
                          "correlation path and were dropped; the band uses %d draws."),
                    B - B_used, B, B_used))

  a <- (1 - level) / 2
  pair_labels <- {
    cp <- utils::combn(K, 2)
    apply(cp, 2, function(ix) sprintf("rho[%d,%d]", ix[1], ix[2]))
  }
  out <- vector("list", C); names(out) <- pair_labels
  for (cc in seq_len(C)) {
    mat <- vapply(paths, function(pc) pc[cc, ], numeric(Tn))  # T x B_used
    lo  <- apply(mat, 1, stats::quantile, probs = a,     names = FALSE)
    hi  <- apply(mat, 1, stats::quantile, probs = 1 - a, names = FALSE)
    out[[cc]] <- cbind(fit = point[cc, ], lower = lo, upper = hi)
  }
  attr(out, "level")  <- level
  attr(out, "B_used") <- B_used
  out
}
