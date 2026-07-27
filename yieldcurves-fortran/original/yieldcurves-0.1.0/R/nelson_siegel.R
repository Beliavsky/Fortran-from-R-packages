#' Fit Nelson-Siegel Yield Curve
#'
#' Estimate a Nelson-Siegel (1987) yield curve model from observed
#' maturity-rate pairs. The model decomposes the yield curve into three
#' factors: level, slope, and curvature.
#'
#' The Nelson-Siegel model is:
#' \deqn{r(m) = \beta_0 + \beta_1 \frac{1 - e^{-m/\tau}}{m/\tau} +
#' \beta_2 \left(\frac{1 - e^{-m/\tau}}{m/\tau} - e^{-m/\tau}\right)}
#'
#' @param maturities Numeric vector of maturities in years.
#' @param rates Numeric vector of observed yields as decimals.
#' @param tau_init Numeric. Initial value for the decay parameter tau.
#'   Default is 1.
#' @param weights Optional numeric vector of weights for each observation.
#'   Must be the same length as `maturities`. Useful for emphasising
#'   liquid tenors. If NULL (default), all observations are equally weighted.
#' @param type Character. Rate type: `"zero"` (default), `"par"`, or
#'   `"forward"`.
#' @param date Optional Date for the curve.
#'
#' @return A `yc_curve` object with `method = "nelson_siegel"` and
#'   `params` containing `beta0`, `beta1`, `beta2`, and `tau`.
#'
#' @references
#' Nelson, C.R. and Siegel, A.F. (1987). Parsimonious Modeling of Yield
#' Curves. *The Journal of Business*, 60(4), 473--489.
#' \doi{10.1086/296409}
#'
#' @export
#' @examples
#' maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
#' rates <- c(0.052, 0.050, 0.048, 0.045, 0.043, 0.042, 0.041,
#'            0.040, 0.042, 0.043)
#' fit <- yc_nelson_siegel(maturities, rates)
#' fit
yc_nelson_siegel <- function(maturities, rates, tau_init = 1,
                              weights = NULL,
                              type = c("zero", "par", "forward"),
                              date = NULL) {
  type <- match.arg(type)
  validate_maturities(maturities)
  validate_rates(rates, length(maturities))
  validate_positive_scalar(tau_init, "tau_init")

  if (!is.null(weights)) {
    if (!is.numeric(weights) || length(weights) != length(maturities)) {
      cli_abort("{.arg weights} must be a numeric vector of length {length(maturities)}.")
    }
    if (any(weights < 0)) {
      cli_abort("{.arg weights} must be non-negative.")
    }
  }

  if (length(maturities) < 4) {
    cli_warn("Nelson-Siegel model has 4 parameters; fitting with fewer than 4 data points may be unreliable.")
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

  # Multi-start grid search over tau (BIS/ECB methodology)
  # For each tau, the NS model is linear in (beta0, beta1, beta2),
  # so solve via weighted OLS to find the best starting point.
  tau_grid <- c(0.5, 1, 2, 3, 5, 8, 12)
  best_sse <- Inf
  best_start <- NULL

  W <- diag(w)

  for (tau_try in tau_grid) {
    X <- ns_loadings(maturities, tau_try)
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
      best_start <- c(beta_ols, tau_try)
    }
  }

  # Fallback if grid search fails
 if (is.null(best_start)) {
    beta0_init <- rates[length(rates)]
    beta1_init <- rates[1] - rates[length(rates)]
    beta2_init <- 2 * mean(rates) - rates[1] - rates[length(rates)]
    best_start <- c(beta0_init, beta1_init, beta2_init, tau_init)
  }

  # Objective function
  obj <- function(par) {
    beta0 <- par[1]
    beta1 <- par[2]
    beta2 <- par[3]
    tau   <- par[4]
    fitted <- ns_rate(maturities, beta0, beta1, beta2, tau)
    sum(w * (rates - fitted)^2)
  }

  # Local optimisation from best grid point
  result <- optim(
    par = best_start,
    fn = obj,
    method = "L-BFGS-B",
    lower = c(-Inf, -Inf, -Inf, 0.01),
    upper = c(Inf, Inf, Inf, 30)
  )

  if (result$convergence != 0) {
    cli_warn("Nelson-Siegel optimisation did not fully converge (code {result$convergence}). Results may be approximate.")
  }

  beta0 <- result$par[1]
  beta1 <- result$par[2]
  beta2 <- result$par[3]
  tau   <- result$par[4]

  fitted_rates <- ns_rate(maturities, beta0, beta1, beta2, tau)

  new_yc_curve(
    maturities = maturities,
    rates = rates,
    type = type,
    method = "nelson_siegel",
    params = list(beta0 = beta0, beta1 = beta1, beta2 = beta2, tau = tau),
    fitted = fitted_rates,
    residuals = rates - fitted_rates,
    date = date
  )
}
