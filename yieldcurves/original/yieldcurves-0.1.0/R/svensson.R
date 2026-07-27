#' Fit Svensson Yield Curve
#'
#' Estimate a Svensson (1994) yield curve model from observed maturity-rate
#' pairs. Extends Nelson-Siegel by adding a second curvature term with its
#' own decay parameter, providing greater flexibility for curves with two
#' humps.
#'
#' @param maturities Numeric vector of maturities in years.
#' @param rates Numeric vector of observed yields as decimals.
#' @param tau1_init Numeric. Initial value for the first decay parameter.
#'   Default is 1.
#' @param tau2_init Numeric. Initial value for the second decay parameter.
#'   Default is 5.
#' @param weights Optional numeric vector of weights for each observation.
#'   Must be the same length as `maturities`. If NULL (default), all
#'   observations are equally weighted.
#' @param type Character. Rate type: `"zero"` (default), `"par"`, or
#'   `"forward"`.
#' @param date Optional Date for the curve.
#'
#' @return A `yc_curve` object with `method = "svensson"` and `params`
#'   containing `beta0`, `beta1`, `beta2`, `beta3`, `tau1`, and `tau2`.
#'
#' @references
#' Svensson, L.E.O. (1994). Estimating and Interpreting Forward Interest
#' Rates: Sweden 1992--1994. *NBER Working Paper*, 4871.
#' \doi{10.3386/w4871}
#'
#' @export
#' @examples
#' maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
#' rates <- c(0.052, 0.050, 0.048, 0.045, 0.043, 0.042, 0.041,
#'            0.040, 0.042, 0.043)
#' fit <- yc_svensson(maturities, rates)
#' fit
yc_svensson <- function(maturities, rates, tau1_init = 1, tau2_init = 5,
                         weights = NULL,
                         type = c("zero", "par", "forward"),
                         date = NULL) {
  type <- match.arg(type)
  validate_maturities(maturities)
  validate_rates(rates, length(maturities))
  validate_positive_scalar(tau1_init, "tau1_init")
  validate_positive_scalar(tau2_init, "tau2_init")

  if (!is.null(weights)) {
    if (!is.numeric(weights) || length(weights) != length(maturities)) {
      cli_abort("{.arg weights} must be a numeric vector of length {length(maturities)}.")
    }
    if (any(weights < 0)) {
      cli_abort("{.arg weights} must be non-negative.")
    }
  }

  if (length(maturities) < 6) {
    cli_warn("Svensson model has 6 parameters; fitting with fewer than 6 data points may be unreliable.")
  }

  if (!is.null(date) && !inherits(date, "Date")) {
    cli_abort("{.arg date} must be a {.cls Date} or NULL.")
  }

  # Sort by maturity
  ord <- order(maturities)
  maturities <- maturities[ord]
  rates <- rates[ord]
  if (!is.null(weights)) weights <- weights[ord]

  # Default weights
  w <- if (is.null(weights)) rep(1, length(maturities)) else weights

  # Multi-start grid search over (tau1, tau2) pairs (BIS/ECB methodology)
  # For each (tau1, tau2), the Svensson model is linear in
  # (beta0, beta1, beta2, beta3), so solve via weighted OLS.
  tau1_grid <- c(1, 2, 5, 8)
  tau2_grid <- c(0.5, 1, 3, 5)
  best_sse <- Inf
  best_start <- NULL

  W <- diag(w)

  for (tau1_try in tau1_grid) {
    for (tau2_try in tau2_grid) {
      if (abs(tau1_try - tau2_try) < 0.1) next
      X <- sv_loadings(maturities, tau1_try, tau2_try)
      XtWX <- crossprod(X, W %*% X)
      XtWy <- crossprod(X, W %*% rates)
      beta_ols <- tryCatch(
        as.numeric(solve(XtWX, XtWy)),
        error = function(e) NULL
      )
      if (is.null(beta_ols)) next
      fitted_try <- as.numeric(X %*% beta_ols)
      sse <- sum(w * (rates - fitted_try)^2)
      if (sse < best_sse) {
        best_sse <- sse
        best_start <- c(beta_ols, tau1_try, tau2_try)
      }
    }
  }

  # Fallback: use NS fit as starting point
  if (is.null(best_start)) {
    ns_fit <- yc_nelson_siegel(maturities, rates, tau_init = tau1_init,
                                weights = weights, type = type, date = date)
    best_start <- c(ns_fit$params$beta0, ns_fit$params$beta1,
                    ns_fit$params$beta2, 0, ns_fit$params$tau, tau2_init)
  }

  obj <- function(par) {
    fitted <- sv_rate(maturities, par[1], par[2], par[3], par[4], par[5], par[6])
    sum(w * (rates - fitted)^2)
  }

  result <- optim(
    par = best_start,
    fn = obj,
    method = "L-BFGS-B",
    lower = c(-Inf, -Inf, -Inf, -Inf, 0.01, 0.01),
    upper = c(Inf, Inf, Inf, Inf, 30, 30)
  )

  if (result$convergence != 0) {
    cli_warn("Svensson optimisation did not fully converge (code {result$convergence}). Results may be approximate.")
  }

  beta0 <- result$par[1]
  beta1 <- result$par[2]
  beta2 <- result$par[3]
  beta3 <- result$par[4]
  tau1  <- result$par[5]
  tau2  <- result$par[6]

  fitted_rates <- sv_rate(maturities, beta0, beta1, beta2, beta3, tau1, tau2)

  new_yc_curve(
    maturities = maturities,
    rates = rates,
    type = type,
    method = "svensson",
    params = list(beta0 = beta0, beta1 = beta1, beta2 = beta2,
                  beta3 = beta3, tau1 = tau1, tau2 = tau2),
    fitted = fitted_rates,
    residuals = rates - fitted_rates,
    date = date
  )
}
