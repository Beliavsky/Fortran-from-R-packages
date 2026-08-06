#' F2 Variance-Spanning Test (Slopes Only)
#'
#' Tests the null \eqn{H_0:\ \delta = 0} that adding test assets does not
#' improve the minimum-variance frontier spanned by the benchmarks (variance
#' spanning). The statistic compares frontier-defining quantities of the
#' augmented (benchmark + test) universe to those of the benchmark subset.
#'
#' @param R1 Numeric matrix of benchmark returns, dimension \eqn{T \times K}.
#' @param R2 Numeric matrix of test-asset returns, dimension \eqn{T \times N}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{pval}}{P-value for the \eqn{F}-statistic under the null.}
#'   \item{\code{stat}}{F2 \eqn{F}-statistic.}
#'   \item{\code{H0}}{Null hypothesis description, \code{"delta = 0"}.}
#' }
#'
#' @details
#' Under standard conditions (i.i.d. returns, full-rank covariances), the reference
#' distribution is \eqn{F_{N,\ T-K-N+1}}. Finite-sample feasibility requires
#' \eqn{T-K-N+1 \ge 1}.
#'
#' @references
#' \insertRef{KanZhou2012}{spantest}
#'
#' @examples
#' set.seed(123)
#' R1 <- matrix(rnorm(300), 100, 3)  # benchmarks: T=100, K=3
#' R2 <- matrix(rnorm(200), 100, 2)  # tests:      T=100, N=2
#' out <- span_f2(R1, R2)
#' out$stat; out$pval; out$H0
#'
#' @family Variance Spanning Tests
#'
#' @importFrom stats pf cov
#' @export
span_f2 <- function(R1, R2) {

  R <- cbind(R1, R2)
  K <- ncol(R1)
  N <- ncol(R2)
  t <- nrow(R)

  df1 <- N
  df2 <- t - K - N + 1
  if (df2 < 1) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "delta = 0"))
  }

  mu  <- matrix(colMeans(R), ncol = 1)
  V   <- cov(R)
  iV  <- tryCatch(solve(V), error = function(e) NULL)
  one <- matrix(1, K + N, 1)

  mu1  <- matrix(colMeans(R1), ncol = 1)
  V1   <- cov(R1)
  iV1  <- tryCatch(solve(V1), error = function(e) NULL)
  one1 <- matrix(1, K, 1)

  if (is.null(iV) || is.null(iV1)) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "delta = 0"))
  }

  a  <- t(mu) %*% iV %*% mu
  b  <- t(mu) %*% iV %*% one
  c_ <- t(one) %*% iV %*% one
  d  <- a * c_ - b^2

  a1 <- t(mu1) %*% iV1 %*% mu1
  b1 <- t(mu1) %*% iV1 %*% one1
  c1 <- t(one1) %*% iV1 %*% one1
  d1 <- a1 * c1 - b1^2

  # Compute F2 statistic
  stat <- (df2 / df1) * (((c_ + d) / (c1 + d1)) * ((1 + a1) / (1 + a)) - 1)
  pval <- pf(stat, df1, df2, lower.tail = FALSE)

  return(list(pval = as.numeric(pval), stat = as.numeric(stat), H0 = "delta = 0"))
}
