#' Britten–Jones Tangency-Portfolio Spanning Test (1999)
#'
#' Tests whether the tangency (maximum Sharpe) portfolio of the augmented
#' universe (benchmarks + test assets) is spanned by the benchmark assets
#' alone. Following Britten–Jones (1999), the statistic arises from a
#' regression of a constant on the raw asset returns (without an intercept)
#' and yields an \eqn{F} test of the tangency-spanning restriction.
#'
#' @param R1 Numeric matrix of benchmark returns, dimension \eqn{T \times K}.
#' @param R2 Numeric matrix of test-asset returns, dimension \eqn{T \times N}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{pval}}{P-value for the \eqn{F}-statistic under the null.}
#'   \item{\code{stat}}{Britten–Jones \eqn{F}-statistic.}
#'   \item{\code{H0}}{Null hypothesis description, \code{"alpha = 0"} (the test
#'   assets do not expand the tangency portfolio spanned by the benchmarks).}
#' }
#'
#' @details
#' Following Britten–Jones (1999), the (unnormalized) tangency-portfolio weights
#' are obtained by regressing a vector of ones on the raw asset returns without an
#' intercept. The null that the benchmarks span the tangency portfolio is the
#' restriction that the test assets carry zero weight, tested by comparing the
#' residual sums of squares of the full and benchmark-only regressions. The
#' statistic has an \eqn{F_{N,\ T-K-N}} reference distribution; finite-sample
#' feasibility requires \eqn{T-K-N \ge 1}.
#'
#' \strong{Equivalence with [span_grs()].} Although derived from a regression of
#' ones on raw returns rather than from intercepts, this statistic is
#' numerically identical to the Gibbons--Ross--Shanken statistic returned by
#' [span_grs()], to machine precision and for every \eqn{(T, K, N)}. The two are
#' the same test of \eqn{\alpha = 0} computed by different routes; agreement
#' between them validates the implementation but is not independent evidence.
#' [span_f1()] is close but genuinely distinct.
#'
#' @references
#' \insertRef{BrittenJones1999}{spantest} \cr
#'
#' @examples
#' set.seed(321)
#' R1 <- matrix(rnorm(300), 100, 3)  # benchmarks: T=100, K=3
#' R2 <- matrix(rnorm(200), 100, 2)  # tests:      T=100, N=2
#' out <- span_bj(R1, R2)
#' out$stat; out$pval; out$H0
#'
#' @family Alpha Spanning Tests
#'
#' @importFrom stats pf
#' @export
span_bj <- function(R1, R2) {
  R <- cbind(R1, R2)     # Combine benchmark (R1) and test (R2)
  t <- nrow(R)
  K <- ncol(R1)          # K = benchmark assets
  N <- ncol(R2)          # N = test assets

  # DF restriction
  if ((t - K - N) < 1) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0"))
  }

  ones <- rep(1, t)

  # Unrestricted regression of a constant on ALL asset returns (no intercept):
  # the OLS coefficients are proportional to the combined tangency-portfolio
  # weights (Britten-Jones, 1999).
  XtX_inv <- tryCatch(solve(crossprod(R)), error = function(e) NULL)
  if (is.null(XtX_inv)) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0"))
  }
  bU   <- XtX_inv %*% crossprod(R, ones)
  rssU <- sum((ones - R %*% bU)^2)

  # Restricted regression: only the benchmark assets carry weight
  # (test-asset weights constrained to zero).
  XtX1_inv <- tryCatch(solve(crossprod(R1)), error = function(e) NULL)
  if (is.null(XtX1_inv)) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0"))
  }
  bR   <- XtX1_inv %*% crossprod(R1, ones)
  rssR <- sum((ones - R1 %*% bR)^2)

  stat <- ((rssR - rssU) / N) / (rssU / (t - K - N))
  pval <- pf(stat, N, t - K - N, lower.tail = FALSE)

  list(pval = as.numeric(pval), stat = as.numeric(stat), H0 = "alpha = 0")
}

