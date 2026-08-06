#' Gungor–Luger Alpha-Only Spanning Test (2016)
#'
#' Tests the null \eqn{H_0:\ \alpha = 0} that benchmark assets span the mean
#' (intercepts) of the test assets. Following Gungor & Luger (2016), the
#' procedure uses a Monte Carlo (MC) test based on an \eqn{F_{\max}} statistic
#' with residual sign-flip simulations, yielding Least-Favorable (LMC) and
#' Balanced (BMC) MC p-values and a three-way decision rule.
#'
#' @param R1 Numeric matrix of benchmark returns, dimension \eqn{T \times K}.
#' @param R2 Numeric matrix of test-asset returns, dimension \eqn{T \times N}.
#' @param control List of options:
#' \describe{
#'   \item{\code{totsim}}{Number of MC simulations (default \code{500}).}
#'   \item{\code{pval_thresh}}{Significance level for decisions (default \code{0.05}).}
#'   \item{\code{do_trace}}{Logical; print progress (default \code{FALSE}).}
#' }
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{pval_LMC}}{Least-Favorable MC p-value.}
#'   \item{\code{pval_BMC}}{Balanced MC p-value.}
#'   \item{\code{stat}}{Observed \eqn{F_{\max}} statistic.}
#'   \item{\code{Decisions}}{Decision code: \code{1} = Accept, \code{0} = Reject, \code{NA} = Inconclusive.}
#'   \item{\code{Decisions_string}}{Text label: \code{"Accept"}, \code{"Reject"}, or \code{"Inconclusive"}.}
#'   \item{\code{H0}}{Null hypothesis description, \code{"alpha = 0"}.}
#' }
#'
#' @details
#' Accept if \code{pval_LMC > alpha}; Reject if \code{pval_BMC <= alpha};
#' otherwise Inconclusive. The subseries sign-flip MC approach is robust to
#' heteroskedasticity, serial dependence, and heavy tails, making it suitable
#' for high-dimensional settings where classical alpha tests (e.g., GRS) may
#' suffer from size distortions.
#'
#' @references
#' \insertRef{GungorLuger2016}{spantest} \cr
#'
#' @examples
#' set.seed(1234)
#' R1 <- matrix(rnorm(300), 100, 3)
#' R2 <- matrix(rnorm(200), 100, 2)
#' out <- span_gl_a(R1, R2, control = list(totsim = 100, do_trace = FALSE))
#' out$Decisions_string; out$pval_LMC; out$pval_BMC
#'
#' @family Alpha Spanning Tests
#'
#' @importFrom stats rnorm runif
#' @export
span_gl_a <- function(R1, R2, control = list()) {

  con <- list(totsim = 500, do_trace = FALSE, pval_thresh = 0.05)
  con[names(control)] <- control
  thresh <- con$pval_thresh

  X <- R1
  Y <- R2

  K <- ncol(X)
  N <- ncol(Y)
  TT <- nrow(X)
  totsim <- con$totsim
  ones <- matrix(1, TT, 1)

  XX <- cbind(ones, X)

  XX_crossprod <- crossprod(XX)
  # Collinear benchmarks make the design singular. The other tests in the
  # package return NA rather than raising in that case (see span_grs); match it.
  Xtemp <- tryCatch(solve(XX_crossprod), error = function(e) NULL)
  if (is.null(Xtemp))
    return(list(pval_LMC = NA_real_, pval_BMC = NA_real_, stat = NA_real_,
                Decisions = NA, Decisions_string = NA_character_,
                H0 = "alpha = 0"))
  Bhat1 <- Xtemp %*% crossprod(XX, Y)

  Ehat1 <- Y - XX %*% Bhat1
  SSRu <- crossprod(Ehat1)

  H <- matrix(0, 1, K + 1)
  H[1, 1] <- 1
  C <- matrix(0, 1, N)

  HXt <- H %*% Xtemp
  HXtHt <- HXt %*% t(H)
  HXtHt_inv <- solve(HXtHt)
  HB_minus_C <- H %*% Bhat1 - C

  Bhat0 <- Bhat1 - Xtemp %*% t(H) %*% HXtHt_inv %*% HB_minus_C
  Ehat0 <- Y - XX %*% Bhat0
  SSRr <- crossprod(Ehat0)

  diag_SSRr <- diag(SSRr)
  diag_SSRu <- diag(SSRu)
  temp <- (diag_SSRr - diag_SSRu) / diag_SSRu
  Fmax_actual <- max(temp)

  LMCstats <- numeric(totsim)
  BMCstats <- numeric(totsim)
  LMCstats[totsim] <- BMCstats[totsim] <- Fmax_actual

  premult <- Xtemp %*% t(H) %*% HXtHt_inv
  Xtemp_XX <- Xtemp %*% t(XX)

  # Sign-flip simulations. The random signs are drawn in R (so the RNG stream is
  # unchanged); the per-simulation restricted/unrestricted SSR and F-max are
  # computed in C++ (gl_sim_stats), streaming one simulation at a time to avoid
  # the large T x (N*nsim) intermediates. The balanced-MC restricted SSR is
  # constant across sign-flips (esim^2 == Ehat0^2), so it is passed once.
  sim_count <- totsim - 1
  sign_mat <- matrix(sign(rnorm(TT * sim_count)), TT, sim_count)
  sim <- gl_sim_stats(XX, Xtemp_XX, XX %*% Bhat0, Ehat0,
                      H, C, premult, sign_mat, colSums(Ehat0^2))
  LMCstats[1:sim_count] <- sim$LMC
  BMCstats[1:sim_count] <- sim$BMC

  uu <- runif(totsim)
  GL_pval_LMC <- (totsim - f_ranklex(LMCstats, uu) + 1) / totsim
  GL_pval_BMC <- (totsim - f_ranklex(BMCstats, uu) + 1) / totsim

  Decisions <- if (GL_pval_LMC > thresh) {
    1
  } else if (GL_pval_BMC <= thresh) {
    0
  } else {
    NA
  }

  Decisions_string <- if (GL_pval_LMC > thresh) {
    "Accept"
  } else if (GL_pval_BMC <= thresh) {
    "Reject"
  } else {
    "Inconclusive"
  }

  if (con$do_trace) {
    cat("============================================\n")
    cat("F-max:", Fmax_actual, "\n")
    cat("LMC p-value:", GL_pval_LMC, "\n")
    cat("BMC p-value:", GL_pval_BMC, "\n")
    cat("Decision:", Decisions_string, "\n")
  }

  GL <- list(
    pval_LMC = GL_pval_LMC,
    pval_BMC = GL_pval_BMC,
    stat = Fmax_actual,
    Decisions = Decisions,
    Decisions_string = Decisions_string,
    H0 = "alpha = 0"
  )

  return(GL)
}
