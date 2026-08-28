#' Calculate EFA Fit Statistics
#'
#' @param x An object of class "GPArotation".
#' @param R Optional. The observed correlation matrix. If NULL extracted from
#'   \code{x$correlation}, which is stored automatically when a \code{factanal}
#'   object is passed as input.
#' @param n.obs Optional. Sample size. If NULL extracted from \code{x$n.obs}.
#'   Required for RMSEA.
#' @param alpha Significance level for RMSEA confidence interval. Default 0.10
#'   giving a 90\% confidence interval, the standard for RMSEA.
#'
#' @return A list containing dof, chisq, SRMR, RMSEA, RMSEA.l, RMSEA.u, alpha.
#' @export
calc_fitstats <- function(x, R = NULL, n.obs = NULL, alpha = 0.10) {

  L   <- unclass(x$loadings)
  Phi <- x$Phi
  v   <- nrow(L)
  k   <- ncol(L)

  # 1. Observed correlation matrix
  if (is.null(R))
    R <- x$correlation

  if (is.null(R))
    stop("Correlation matrix not found. Pass a factanal object to store it ",
         "automatically, or supply R manually.")

  # 2. Sample size
  if (is.null(n.obs))
    n.obs <- x$n.obs

  # 3. Degrees of freedom
  dof <- ((v - k)^2 - (v + k)) / 2

  if (dof <= 0)
    warning("Model is saturated or underidentified (df <= 0). ",
            "Fit statistics may be invalid.")

  # 4. Model-implied correlation matrix
  if (is.null(Phi))
    Phi <- diag(k)

  R_hat       <- L %*% Phi %*% t(L)
  diag(R_hat) <- 1

  # 5. Residual matrix and SRMR
  Res       <- R - R_hat
  lower_res <- Res[lower.tri(Res, diag = FALSE)]
  srmr      <- sqrt(mean(lower_res^2))

  # 6. ML fit function and chi-square
  # F_min = log|R_hat| + tr(R R_hat^-1) - log|R| - v
  # Strictly correct for MLE extraction; approximation for PAF/minres
  F_min     <- log(det(R_hat)) + sum(diag(R %*% solve(R_hat))) -
               log(det(R)) - v
  chisq_emp <- F_min * (n.obs - 1)

  # 7. RMSEA - matching psych formula
  rmsea   <- NA
  rmsea_l <- NA
  rmsea_u <- NA

  if (!is.null(n.obs) && dof > 0) {
    rmsea <- sqrt(max(chisq_emp / (dof * n.obs) - 1 / (n.obs - 1), 0))

    tail    <- alpha / 2
    max_ncp <- max(n.obs, chisq_emp) + 2 * n.obs

    # Lower bound
    rmsea_l <- 0
    if (pchisq(chisq_emp, df = dof) > (1 - tail)) {
      rmsea_l <- tryCatch(
        sqrt(uniroot(function(ncp)
          pchisq(chisq_emp, df = dof, ncp = ncp) - (1 - tail),
          c(0, max_ncp))$root / ((n.obs - 1) * dof)),
        error = function(e) NA
      )
    }

    # Upper bound
    rmsea_u <- 0
    if (pchisq(chisq_emp, df = dof) > tail) {
      rmsea_u <- tryCatch(
        sqrt(uniroot(function(ncp)
          pchisq(chisq_emp, df = dof, ncp = ncp) - tail,
          c(0, max_ncp))$root / ((n.obs - 1) * dof)),
        error = function(e) {
          if (rmsea <= 0) 0 else {
            message("Could not find RMSEA upper bound.")
            NA
          }
        }
      )
    }
  }

  list(
    dof     = dof,
    chisq   = chisq_emp,
    SRMR    = srmr,
    RMSEA   = rmsea,
    RMSEA.l = rmsea_l,
    RMSEA.u = rmsea_u,
    alpha   = alpha
  )
}