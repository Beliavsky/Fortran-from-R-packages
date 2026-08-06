#' Kempf–Memmel GMVP Spanning Test
#'
#' Tests whether the Global Minimum Variance Portfolio (GMVP) of the combined
#' (benchmark + test) universe equals the GMVP of the benchmark assets alone.
#' Following Kempf & Memmel (2006), the null assesses whether adding new assets
#' improves the minimum-variance frontier.
#'
#' @param R1 Numeric matrix of benchmark returns, dimension \eqn{T \times K}.
#' @param R2 Numeric matrix of test-asset returns, dimension \eqn{T \times N}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{pval}}{P-value for the F-statistic under the null.}
#'   \item{\code{stat}}{F-statistic value.}
#'   \item{\code{H0}}{Null hypothesis description, \code{"delta = 0"} (equivalently,
#'   GMVP of the benchmark set equals GMVP of the full universe).}
#' }
#'
#' @details
#' The null hypothesis \eqn{H_0} is that augmenting the benchmark set with the
#' test assets does not change the GMVP weights (\eqn{\Delta = 0}), i.e.,
#' the GMVP of the full universe coincides with that of the benchmark subset.
#' The test uses the Kempf--Memmel regression
#' \deqn{r_{1t} = \alpha + \sum_{j>1} \beta_j (r_{1t} - r_{jt}) + \varepsilon_t,}
#' whose slopes are the global-minimum-variance portfolio weights
#' (\eqn{w_j = \beta_j} for \eqn{j>1} and \eqn{w_1 = 1 - \sum_{j>1}\beta_j}).
#' The test assets leave the GMVP unimproved exactly when their \eqn{\beta}
#' block is zero, so \eqn{H_0} is an \eqn{N}-dimensional linear restriction,
#' tested by an \eqn{F} statistic with reference distribution
#' \eqn{F_{N,\ T-K-N}}. Finite-sample feasibility requires \eqn{T-K-N \ge 1}.
#'
#' Note that regressing a vector of ones on the returns instead---the
#' Britten--Jones tangency construction used by [span_bj()]---produces a test of
#' \eqn{\alpha = 0}, not of \eqn{\Delta = 0}. The two are easily confused: they
#' share the difference-of-returns design matrix and differ only in the
#' dependent variable and the intercept.
#'
#' @references
#' \insertRef{KempfMemmel2006}{spantest} \cr
#'
#' @examples
#' set.seed(123)
#' R1 <- matrix(rnorm(300), 100, 3)  # benchmarks: T=100, K=3
#' R2 <- matrix(rnorm(200), 100, 2)  # test assets: T=100, N=2
#' ans <- span_km(R1, R2)
#' ans$pval; ans$stat; ans$H0
#'
#' @family Variance Spanning Tests
#'
#' @importFrom stats pf
#' @export
span_km <- function(R1, R2) {

  R <- cbind(R1, R2)
  TT <- nrow(R)
  p <- ncol(R1)
  p2 <- ncol(R2)

  # DF restriction
  if ((TT - p - p2) < 1) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "delta = 0"))
  }

  # Kempf-Memmel (2006) GMVP regression:
  #     r_1t = alpha + sum_{j>1} beta_j (r_1t - r_jt) + eps_t,
  # whose slopes ARE the global-minimum-variance weights (w_j = beta_j for
  # j > 1, w_1 = 1 - sum beta_j). The test assets leave the GMVP unimproved
  # exactly when their beta block is zero, which is the delta = 0 null.
  #
  # NOTE. Regressing a vector of ones on the returns instead -- the
  # Britten-Jones tangency construction -- yields an alpha test, not this one.
  # An earlier implementation did that: it had 100% power against alpha and
  # none at all against delta. The regression below is the delta test.
  y <- R[, 1]

  Diff <- sweep(R[, -1, drop = FALSE], 1, R[, 1], FUN = function(x, y) y - x)

  X <- cbind(1, Diff)                       # intercept + (K + N - 1) differences
  XtX <- crossprod(X)
  XtX_inv <- tryCatch(solve(XtX), error = function(e) return(NULL))
  if (is.null(XtX_inv)) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "delta = 0"))
  }

  beta_hat <- XtX_inv %*% crossprod(X, y)
  residuals <- y - X %*% beta_hat
  sigma2 <- drop(crossprod(residuals) / (TT - ncol(X)))

  # Select the test-asset slope block: X is [intercept | (p-1) benchmark
  # differences | p2 test-asset differences], so the block starts after p columns.
  offset <- p
  C <- cbind(matrix(0, p2, offset), diag(p2))

  middle <- tryCatch(solve(C %*% XtX_inv %*% t(C)), error = function(e) return(NULL))
  if (is.null(middle)) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "delta = 0"))
  }

  theta_part <- C %*% beta_hat
  F_stat <- as.numeric(t(theta_part) %*% middle %*% theta_part / (p2 * sigma2))
  p_val <- pf(F_stat, p2, TT - ncol(X), lower.tail = FALSE)

  return(list(pval = p_val, stat = F_stat, H0 = "delta = 0"))
}

