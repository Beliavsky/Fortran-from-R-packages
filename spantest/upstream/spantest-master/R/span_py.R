#' Pesaran–Yamagata Alpha Spanning Test (2024)
#'
#' Implements the Pesaran–Yamagata test for the joint null that all intercepts
#' are zero in a multi-factor spanning regression with possible cross-sectional
#' dependence across test assets.
#'
#' @param R1 Numeric matrix of benchmark returns, dimension \eqn{T \times K}.
#' @param R2 Numeric matrix of test-asset returns, dimension \eqn{T \times N}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{pval}}{P-value under the standard normal reference distribution.}
#'   \item{\code{stat}}{Standardized test statistic.}
#'   \item{\code{H0}}{Null hypothesis description, \code{"alpha = 0"}.}
#' }
#'
#' @details
#' The null hypothesis is that all intercepts are zero (\eqn{\alpha = 0}), meaning
#' the benchmark assets span the expected returns of the test assets. The statistic
#' adjusts for cross-sectional dependence via the residual covariance and has an
#' asymptotic \eqn{\mathcal{N}(0,1)} reference under large \eqn{T,N}. Finite-sample
#' safeguards require \eqn{T-K-1 > 4} and at least two test assets (\eqn{N \ge 2});
#' otherwise \code{pval} and \code{stat} are returned as \code{NA}.
#'
#' @references
#' \insertRef{PesaranYamagata2024}{spantest}
#'
#' @examples
#' set.seed(123)
#' R1 <- matrix(rnorm(300), 100, 3)  # benchmarks: T=100, K=3
#' R2 <- matrix(rnorm(200), 100, 2)  # tests:      T=100, N=2
#' out <- span_py(R1, R2)
#' out$pval; out$stat; out$H0
#'
#' @family Alpha Spanning Tests
#'
#' @importFrom stats pnorm qnorm
#' @export
span_py <- function(R1, R2) {
  X <- R1
  Y <- R2
  t <- nrow(X)
  N <- ncol(Y)
  K <- ncol(X)

  # The PY statistic averages pairwise residual correlations, so it requires
  # at least two test assets (N >= 2); with N = 1 the pairwise sum is empty
  # and 0.05 / (N - 1) divides by zero.
  if (N < 2 || (t - K - N) < 1 || (t - K - 1) <= 4) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0"))
  }

  XX <- cbind(1, X)
  XtX_inv <- tryCatch(solve(crossprod(XX)), error = function(e) NULL)
  X_crossprod_inv <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(XtX_inv) || is.null(X_crossprod_inv)) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0"))
  }

  Bhat1 <- XtX_inv %*% crossprod(XX, Y)
  Ehat1 <- Y - XX %*% Bhat1
  SigmaU <- crossprod(Ehat1) / t

  v <- t - K - 1
  ones <- matrix(1, t, 1)
  X_ones <- X %*% X_crossprod_inv %*% crossprod(X, ones)
  MX_ones <- ones - X_ones

  num_scalar <- crossprod(ones, MX_ones)[1, 1] * v
  num <- matrix(num_scalar, N, 1)
  diag_SigmaU <- diag(SigmaU)
  t2 <- (Bhat1[1, ]^2) * num / (t * diag_SigmaU)

  pN <- 0.05 / (N - 1)
  thetaN <- qnorm(1 - pN / 2)^2
  rhobar <- 0
  for (i in 2:N) {
    for (j in 1:(i - 1)) {
      temp <- SigmaU[i, j] / sqrt(SigmaU[i, i] * SigmaU[j, j])
      temp2 <- temp^2
      if (v * temp2 >= thetaN) {
        rhobar <- rhobar + temp2
      }
    }
  }
  rhobar <- rhobar * 2 / (N * (N - 1))

  Jalpha2 <- sum(t2 - v / (v - 2)) / sqrt(N)
  den <- (v / (v - 2)) * sqrt(2 * (v - 1) * (1 + (N - 1) * rhobar) / (v - 4))
  stat <- Jalpha2 / den
  pval <- pnorm(stat, lower.tail = FALSE)

  list(pval = as.numeric(pval), stat = as.numeric(stat), H0 = "alpha = 0")
}

