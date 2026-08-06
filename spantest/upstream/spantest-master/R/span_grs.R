#' Gibbons–Ross–Shanken (GRS) Alpha Spanning Test (1989)
#'
#' Implements the GRS test of the joint null \eqn{H_0:\ \alpha = 0} in the
#' multivariate regression of test-asset returns on benchmark portfolios
#' (with an intercept). The statistic assumes homoskedastic, normally
#' distributed errors and is most reliable when \eqn{T} is large relative
#' to \eqn{K} (benchmarks) and \eqn{N} (test assets).
#'
#' @param R1 Numeric matrix of benchmark returns, dimension \eqn{T \times K}.
#' @param R2 Numeric matrix of test-asset returns, dimension \eqn{T \times N}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{pval}}{P-value for the \eqn{F}-statistic under the null.}
#'   \item{\code{stat}}{GRS \eqn{F}-statistic.}
#'   \item{\code{H0}}{Null hypothesis description, \code{"alpha = 0"}.}
#' }
#'
#' @details
#' Under standard conditions, the reference distribution is
#' \eqn{F_{N,\ T-N-K}}. Finite-sample feasibility requires \eqn{T-N-K \ge 1}.
#'
#' The residual covariance enters in its maximum-likelihood form (divisor
#' \eqn{T}), which is what pairs with the \eqn{(T-N-K)/N} scaling to give an
#' exact \eqn{F} test.
#'
#' \strong{Equivalence with [span_bj()].} This statistic is numerically
#' identical to the Britten--Jones (1999) statistic returned by [span_bj()],
#' to machine precision and for every \eqn{(T, K, N)}: the two are different
#' computational routes to the same test of \eqn{\alpha = 0}, not independent
#' procedures. Reporting both is therefore a check on the implementation, not
#' corroborating evidence from two tests. [span_f1()] is close but genuinely
#' distinct.
#'
#' @references
#' \insertRef{GRS1989}{spantest} \cr
#'
#' @examples
#' set.seed(42)
#' R1 <- matrix(rnorm(300), 100, 3)  # benchmarks: T=100, K=3
#' R2 <- matrix(rnorm(200), 100, 2)  # tests:      T=100, N=2
#' out <- span_grs(R1, R2)
#' out$stat; out$pval; out$H0
#'
#' @family Alpha Spanning Tests
#'
#' @importFrom stats pf
#' @export
span_grs <- function(R1, R2) {

  X <- R1  # Benchmark assets
  Y <- R2  # Test assets

  t <- nrow(X)
  K <- ncol(X)
  N <- ncol(Y)

  if ((t - N - K) < 1) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0"))
  }

  ones <- matrix(1, t, 1)
  XX <- cbind(ones, X)

  # Return NA (rather than error) if any inverse is singular, as the other tests do.
  stat <- tryCatch({
    Bhat1 <- solve(crossprod(XX)) %*% crossprod(XX, Y)
    SSRu <- crossprod(Y - XX %*% Bhat1)
    SigmaU <- SSRu / t

    ahat <- Bhat1[1, ]
    mufactor <- matrix(colMeans(X), K, 1)
    ssqm <- crossprod(scale(X, center = TRUE, scale = FALSE)) / t
    mufactor_ssqm <- solve(ssqm, mufactor)

    denom <- 1 + crossprod(mufactor, mufactor_ssqm)
    num <- ((t - N - K) / N) * crossprod(ahat, solve(SigmaU, ahat))
    as.numeric(num / denom)
  }, error = function(e) NA_real_)

  if (is.na(stat)) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0"))
  }
  pval <- pf(stat, N, t - N - K, lower.tail = FALSE)

  list(pval = pval, stat = stat, H0 = "alpha = 0")
}
