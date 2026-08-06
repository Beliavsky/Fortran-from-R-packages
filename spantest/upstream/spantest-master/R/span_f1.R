#' F1 Alpha-Spanning Test (Intercepts Only)
#'
#' Tests the null \eqn{H_0:\ \alpha = 0} that the intercepts of the test assets
#' are jointly zero when regressed on the benchmark assets, i.e., benchmarks
#' span the mean of the test assets. This is the F1 test of Kan & Zhou (2012).
#'
#' @param R1 Numeric matrix of benchmark returns, dimension \eqn{T \times K}.
#' @param R2 Numeric matrix of test-asset returns, dimension \eqn{T \times N}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{pval}}{P-value for the \eqn{F}-statistic under the null.}
#'   \item{\code{stat}}{F1 \eqn{F}-statistic.}
#'   \item{\code{H0}}{Null hypothesis description, \code{"alpha = 0"}.}
#' }
#'
#' @details
#' Under standard assumptions (i.i.d. returns, full-rank covariances), the
#' reference distribution is \eqn{F_{N,\ T-K-N}}. Finite-sample feasibility
#' requires \eqn{T-K-N \ge 1}.
#'
#' @references
#' \insertRef{KanZhou2012}{spantest}
#'
#' @examples
#' set.seed(123)
#' R1 <- matrix(rnorm(300), 100, 3)  # benchmarks: T=100, K=3
#' R2 <- matrix(rnorm(200), 100, 2)  # tests:      T=100, N=2
#' out <- span_f1(R1, R2)
#' out$stat; out$pval; out$H0
#'
#' @family Alpha Spanning Tests
#'
#' @importFrom stats pf cov
#' @export
span_f1 <- function(R1, R2) {
  R <- cbind(R1, R2)
  TT <- nrow(R)
  K <- ncol(R1)
  N <- ncol(R2)

  df1 <- N
  df2 <- TT - K - N
  if (df2 < 1) return(list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0"))

  mu_full <- matrix(colMeans(R), ncol = 1)
  mu_bench <- matrix(colMeans(R1), ncol = 1)

  Sigma_full <- cov(R)
  Sigma_bench <- cov(R1)

  iSigma_full <- tryCatch(solve(Sigma_full), error = function(e) return(NULL))
  iSigma_bench <- tryCatch(solve(Sigma_bench), error = function(e) return(NULL))
  if (is.null(iSigma_full) || is.null(iSigma_bench)) return(list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0"))

  a <- t(mu_full) %*% iSigma_full %*% mu_full
  a1 <- t(mu_bench) %*% iSigma_bench %*% mu_bench

  stat <- (df2 / df1) * ((a - a1) / (1 + a1))
  pval <- pf(stat, df1, df2, lower.tail = FALSE)

  list(pval = as.numeric(pval), stat = as.numeric(stat), H0 = "alpha = 0")
}
