#' Fit an Extended Vasicek Credit Loss Model
#'
#' Fits the extended Vasicek single-factor credit loss model where the
#' probability of default depends on macroeconomic covariates. All parameters,
#' including asset value correlation, are estimated via closed-form
#' probit-transformed OLS.
#'
#' @param formula a formula of the form \code{response ~ predictors} where the
#'   response is an observed default/loss rate in (0, 1) (or [0, 1] when
#'   \code{portfolio_size} is supplied). The model must include an intercept.
#' @param data a data frame containing the variables in \code{formula}.
#' @param bias_correct logical; if \code{TRUE}, apply the small-sample bias
#'   correction to the variance estimate by multiplying by \eqn{N/(N - m - 1)}.
#' @param portfolio_size optional positive integer. If supplied, a
#'   finite-portfolio variance correction is applied to the response before
#'   fitting (see Yang, 2014, section 4.3). The correction shrinks the
#'   observed rates towards their mean, mapping rates of exactly 0 or 1 into
#'   (0, 1), so such observations are permitted when \code{portfolio_size}
#'   is given.
#'
#' @return An object of class \code{"vasicekfit"}, which is a list containing:
#' \describe{
#'   \item{p}{estimated probability of default (baseline PD)}
#'   \item{rho}{estimated asset value correlation}
#'   \item{kappa}{named numeric vector of macro-factor sensitivities}
#'   \item{sigma2}{MLE variance estimate in probit space}
#'   \item{lm_fit}{the underlying \code{\link{lm}} object}
#'   \item{fitted.values}{fitted values in probit (Y) space}
#'   \item{residuals}{residuals in probit (Y) space}
#'   \item{formula}{the model formula}
#'   \item{call}{the matched call}
#'   \item{terms}{the \code{\link{terms}} object}
#'   \item{model}{the model frame}
#'   \item{bias_correct}{logical flag used}
#'   \item{portfolio_size}{portfolio size used, or \code{NULL}}
#' }
#'
#' @references
#' Vasicek, O. A. (2002). The distribution of loan portfolio value.
#' \emph{Risk}, 15(12), 160--162.
#'
#' Yang, B. H. (2014). Estimating Long-Run PD, Asset Correlation, and
#' Portfolio Level PD by Vasicek Models. MPRA Paper No. 57244.
#'
#' @examples
#' set.seed(42)
#' n <- 100
#' u1 <- rnorm(n)
#' u2 <- rnorm(n)
#' p_true <- 0.03
#' rho_true <- 0.02
#' kappa_true <- c(0.13, -0.07)
#' z <- rnorm(n)
#' x <- stats::pnorm(
#'   (stats::qnorm(p_true) + kappa_true[1] * u1 + kappa_true[2] * u2 +
#'     sqrt(rho_true) * z) / sqrt(1 - rho_true)
#' )
#' d <- data.frame(default_rate = x, unemp = u1, hpi = u2)
#' fit <- vasicekfit(default_rate ~ unemp + hpi, data = d)
#' fit
#'
#' @export
vasicekfit <- function(formula, data, bias_correct = FALSE,
                       portfolio_size = NULL) {

  cl <- match.call()
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  mt <- stats::terms(mf)

  if (attr(mt, "intercept") == 0L) {
    stop("the model must include an intercept; remove '- 1' or '+ 0' from the formula")
  }

  if (!is.numeric(y)) {
    stop("response must be numeric")
  }

  if (!is.null(portfolio_size)) {
    if (portfolio_size < 2) stop("portfolio_size must be >= 2")
    if (any(y < 0 | y > 1)) stop("response must have values in [0, 1]")
    p0 <- mean(y)
    vr <- stats::var(y)
    v0 <- vr - (p0 * (1 - p0) - vr) / (portfolio_size - 1)
    if (v0 <= 0) stop("finite-portfolio correction yielded non-positive variance")
    y <- p0 + (y - p0) * sqrt(v0 / vr)
  }

  if (any(y <= 0 | y >= 1)) {
    stop("response must be numeric with values in (0, 1)",
         if (is.null(portfolio_size))
           "; observed rates of exactly 0 or 1 can be handled by supplying portfolio_size"
         else
           " after the finite-portfolio correction")
  }

  Y <- stats::qnorm(y)

  X <- stats::model.matrix(mt, mf)
  df_ols <- data.frame(Y = Y, X[, -1, drop = FALSE])
  lm_fit <- stats::lm(Y ~ ., data = df_ols)

  beta_hat <- stats::coef(lm_fit)
  N <- length(Y)
  sigma2 <- mean(lm_fit$residuals^2)

  if (bias_correct) {
    m <- ncol(X) - 1L
    sigma2 <- sigma2 * N / (N - m - 1L)
  }

  denom <- sqrt(1 + sigma2)
  p_hat <- stats::pnorm(beta_hat[1L] / denom)
  rho_hat <- sigma2 / (1 + sigma2)

  predictor_names <- names(beta_hat)[-1L]
  kappa_hat <- beta_hat[-1L] / denom
  names(kappa_hat) <- predictor_names

  structure(
    list(
      p              = unname(p_hat),
      rho            = rho_hat,
      kappa          = kappa_hat,
      sigma2         = sigma2,
      lm_fit         = lm_fit,
      fitted.values  = lm_fit$fitted.values,
      residuals      = lm_fit$residuals,
      formula        = formula,
      call           = cl,
      terms          = mt,
      model          = mf,
      bias_correct   = bias_correct,
      portfolio_size = portfolio_size
    ),
    class = "vasicekfit"
  )
}
