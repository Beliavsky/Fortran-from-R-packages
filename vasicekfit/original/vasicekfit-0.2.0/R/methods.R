#' @export
print.vasicekfit <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nVasicek model parameters:\n")
  cat("  p (PD)  =", format(x$p, digits = 6), "\n")
  cat("  rho (AVC) =", format(x$rho, digits = 6), "\n")
  if (length(x$kappa) > 0) {
    cat("  kappa:\n")
    print(x$kappa, digits = 6)
  }
  invisible(x)
}

#' @export
summary.vasicekfit <- function(object, type = c("iid", "HAC"), ...) {
  type <- match.arg(type)

  cat("\nCall:\n")
  print(object$call)

  V <- stats::vcov(object, type = type, ...)
  cf <- stats::coef(object)
  ses <- sqrt(diag(V))
  zval <- cf / ses
  pval <- 2 * stats::pnorm(-abs(zval))

  # H0: param = 0 is unreachable for p and rho (both strictly positive by
  # construction), so the Wald z / p-value are not interpretable there.
  zval[c("p", "rho")] <- NA_real_
  pval[c("p", "rho")] <- NA_real_

  ctab <- cbind(
    Estimate    = cf,
    `Std. Error` = ses,
    `z value`   = zval,
    `Pr(>|z|)`  = pval
  )

  cat("\nCoefficients (original-space parameters):\n")
  stats::printCoefmat(ctab, has.Pvalue = TRUE, na.print = "")

  cat("\nsigma2 =", format(object$sigma2, digits = 6), "\n")
  cat("Bias correction:", object$bias_correct, "\n")
  if (!is.null(object$portfolio_size)) {
    cat("Portfolio size:", object$portfolio_size, "\n")
  }
  cat("Variance type:", type, "\n")

  cat("\n--- Underlying OLS regression (probit space) ---\n")
  print(summary(object$lm_fit))

  invisible(object)
}

#' @export
vcov.vasicekfit <- function(object, type = c("iid", "HAC"), ...) {
  type <- match.arg(type)
  if (type == "HAC" && !requireNamespace("sandwich", quietly = TRUE)) {
    stop("type = \"HAC\" requires the 'sandwich' package; ",
         "install it with install.packages(\"sandwich\").",
         call. = FALSE)
  }

  beta_hat <- stats::coef(object$lm_fit)
  sigma2   <- object$sigma2
  s        <- sqrt(1 + sigma2)
  N        <- length(object$residuals)
  k        <- length(beta_hat)
  m        <- k - 1L

  if (type == "iid") {
    Vbeta <- stats::vcov(object$lm_fit)
    if (isTRUE(object$bias_correct)) {
      var_sigma2 <- 2 * sigma2^2 / (N - m - 1L)
    } else {
      var_sigma2 <- 2 * sigma2^2 * (N - m - 1L) / N^2
    }
    V <- matrix(0, k + 1L, k + 1L)
    V[seq_len(k), seq_len(k)] <- Vbeta
    V[k + 1L, k + 1L] <- var_sigma2
  } else {
    X <- stats::model.matrix(object$lm_fit)
    e <- object$residuals
    A_inv <- solve(crossprod(X) / N)
    ifx_beta <- (X * e) %*% t(A_inv)
    ifx_sig2 <- e^2 - sigma2
    psi <- cbind(ifx_beta, ifx_sig2)
    V <- sandwich::lrvar(psi, ...)
    V <- (V + t(V)) / 2
  }

  beta0  <- beta_hat[1L]
  phi_b0 <- stats::dnorm(beta0 / s)

  J <- matrix(0, k + 1L, k + 1L)
  J[1L, 1L]      <- phi_b0 / s
  J[1L, k + 1L] <- -phi_b0 * beta0 / (2 * s^3)
  J[2L, k + 1L] <- 1 / (1 + sigma2)^2
  if (m > 0L) {
    for (j in seq_len(m)) {
      J[2L + j, 1L + j]  <- 1 / s
      J[2L + j, k + 1L] <- -beta_hat[1L + j] / (2 * s^3)
    }
  }

  out <- J %*% V %*% t(J)
  nms <- c("p", "rho", names(object$kappa))
  dimnames(out) <- list(nms, nms)
  out
}

#' @export
confint.vasicekfit <- function(object, parm, level = 0.95,
                                type = c("iid", "HAC"), ...) {
  type <- match.arg(type)
  cf  <- stats::coef(object)
  ses <- sqrt(diag(stats::vcov(object, type = type, ...)))
  pnms <- names(cf)

  if (missing(parm)) {
    parm <- pnms
  } else if (is.numeric(parm)) {
    parm <- pnms[parm]
  }

  a <- (1 - level) / 2
  z <- stats::qnorm(1 - a)

  ci <- cbind(cf[parm] - z * ses[parm], cf[parm] + z * ses[parm])
  colnames(ci) <- sprintf("%.1f %%", 100 * c(a, 1 - a))
  rownames(ci) <- parm
  ci
}

#' @export
coef.vasicekfit <- function(object, ...) {
  c(p = object$p, rho = object$rho, object$kappa)
}

#' @export
fitted.vasicekfit <- function(object, ...) {
  object$fitted.values
}

#' @export
residuals.vasicekfit <- function(object, ...) {
  object$residuals
}
