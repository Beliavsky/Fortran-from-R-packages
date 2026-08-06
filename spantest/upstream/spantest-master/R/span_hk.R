#' Huberman–Kandel Joint Mean–Variance Spanning Test (1987)
#'
#' Tests the joint null \eqn{H_0:\ \alpha = 0,\ \delta = 0} that the benchmark
#' assets span the mean–variance frontier of the augmented (benchmark + test)
#' universe. Following Huberman & Kandel (1987), the statistic compares the
#' frontiers with and without the additional assets.
#'
#' @param R1 Numeric matrix of benchmark returns, dimension \eqn{T \times K}.
#' @param R2 Numeric matrix of test-asset returns, dimension \eqn{T \times N}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{pval}}{P-value for the \eqn{F}-statistic under the null.}
#'   \item{\code{stat}}{\eqn{F}-statistic value.}
#'   \item{\code{H0}}{Null hypothesis description, \code{"alpha = 0 and delta = 0"}.}
#' }
#'
#' @details
#' The test evaluates whether adding the test assets changes the efficient
#' frontier implied by the benchmarks. Under standard regularity conditions,
#' the statistic has an \eqn{F} reference with \eqn{(2N,\ 2(T-K-N))} degrees of
#' freedom. Finite-sample feasibility requires \eqn{T-K-N \ge 1}.
#'
#' @references
#' \insertRef{HubermanKandel1987}{spantest} \cr
#'
#' @examples
#' set.seed(123)
#' R1 <- matrix(rnorm(300), 100, 3)  # benchmarks: T=100, K=3
#' R2 <- matrix(rnorm(200), 100, 2)  # tests:      T=100, N=2
#' out <- span_hk(R1, R2)
#' out$stat; out$pval; out$H0
#'
#' @family Joint Mean-Variance Spanning Tests
#'
#' @importFrom stats pf cov
#' @export
span_hk <- function(R1, R2) {
  R <- cbind(R1, R2)
  TT <- nrow(R)
  K <- ncol(R1)
  N <- ncol(R2)

  df1 <- 2 * N
  df2 <- 2 * (TT - K - N)
  if (df2 < 1) return(list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0 and delta = 0"))

  mu <- matrix(colMeans(R), ncol = 1)
  mu1 <- matrix(colMeans(R1), ncol = 1)

  Sigma <- cov(R)
  Sigma1 <- cov(R1)
  iSigma <- tryCatch(solve(Sigma), error = function(e) return(NULL))
  iSigma1 <- tryCatch(solve(Sigma1), error = function(e) return(NULL))
  if (is.null(iSigma) || is.null(iSigma1)) return(list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0 and delta = 0"))

  a <- t(mu) %*% iSigma %*% mu
  b <- t(mu) %*% iSigma %*% rep(1, K + N)
  cc <- t(rep(1, K + N)) %*% iSigma %*% rep(1, K + N)
  d <- a * cc - b^2

  a1 <- t(mu1) %*% iSigma1 %*% mu1
  b1 <- t(mu1) %*% iSigma1 %*% rep(1, K)
  c1 <- t(rep(1, K)) %*% iSigma1 %*% rep(1, K)
  d1 <- a1 * c1 - b1^2

  U <- (c1 + d1) / (cc + d)
  stat <- ((TT - K - N) / N) * (1 / sqrt(U) - 1)
  pval <- pf(stat, df1, df2, lower.tail = FALSE)

  list(pval = as.numeric(pval), stat = as.numeric(stat), H0 = "alpha = 0 and delta = 0")
}

