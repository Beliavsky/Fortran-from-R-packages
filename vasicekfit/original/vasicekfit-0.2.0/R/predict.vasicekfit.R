#' Predict Method for vasicekfit Objects
#'
#' Obtain predictions from a fitted Vasicek model, optionally at specified
#' confidence levels.
#'
#' @param object a \code{\link{vasicekfit}} object.
#' @param newdata an optional data frame of new predictor values. If omitted,
#'   the training data are used.
#' @param type character; \code{"link"} (default) returns predictions in probit
#'   space, \code{"response"} returns the conditional mean loss rate (the
#'   effective PD) on the (0, 1) scale.
#' @param alpha optional numeric vector of confidence levels in (0, 1). When
#'   supplied, predictions are conditional quantiles of the loss distribution
#'   at each confidence level. The result is a matrix with one column per
#'   alpha value.
#' @param ... additional arguments (currently unused).
#'
#' @details
#' With \code{type = "response"} and \code{alpha = NULL}, the prediction is
#' the conditional mean of the loss distribution given the covariates (the
#' effective PD), \eqn{\Phi(\Phi^{-1}(p) + \sum_j \kappa_j u_j)}, matching
#' the convention of \code{\link{predict.glm}}. Note that this is not
#' \code{pnorm} of the \code{"link"} prediction: the model has an error term
#' in probit space and the probit transform is nonlinear, so the
#' back-transformed linear predictor \eqn{\Phi(\hat\eta)} is the conditional
#' \emph{median} of the loss distribution (available as \code{alpha = 0.5}),
#' which for \eqn{p < 0.5} is smaller than the mean.
#'
#' @return If \code{alpha} is \code{NULL}, a numeric vector of predictions. If
#'   \code{alpha} is supplied, a matrix with rows corresponding to observations
#'   and columns to confidence levels.
#'
#' @examples
#' set.seed(42)
#' n <- 100
#' u <- rnorm(n)
#' x <- stats::pnorm((stats::qnorm(0.03) + 0.1 * u +
#'   sqrt(0.02) * rnorm(n)) / sqrt(1 - 0.02))
#' d <- data.frame(default_rate = x, macro = u)
#' fit <- vasicekfit(default_rate ~ macro, data = d)
#'
#' # Conditional mean PD (expected loss rate given the covariates)
#' predict(fit, type = "response")
#'
#' # Conditional median loss rate
#' predict(fit, alpha = 0.5, type = "response")
#'
#' # 99th percentile loss rate under stress
#' predict(fit, newdata = data.frame(macro = 2), alpha = 0.99,
#'         type = "response")
#'
#' @export
predict.vasicekfit <- function(object, newdata = NULL,
                               type = c("link", "response"),
                               alpha = NULL, ...) {

  type <- match.arg(type)

  if (is.null(alpha) && type == "link") {
    if (is.null(newdata)) {
      return(object$fitted.values)
    }
    tt <- stats::delete.response(object$terms)
    mf <- stats::model.frame(tt, newdata)
    X <- stats::model.matrix(tt, mf)
    beta <- stats::coef(object$lm_fit)
    return(drop(X %*% beta))
  }

  if (!is.null(alpha) && any(alpha <= 0 | alpha >= 1)) {
    stop("alpha must be in (0, 1)")
  }

  p   <- object$p
  rho <- object$rho
  kappa <- object$kappa

  if (is.null(newdata)) {
    mf_pred <- object$model
  } else {
    tt <- stats::delete.response(object$terms)
    mf_pred <- stats::model.frame(tt, newdata)
  }
  X_pred <- stats::model.matrix(stats::delete.response(object$terms), mf_pred)
  u_contrib <- drop(X_pred[, -1, drop = FALSE] %*% kappa)

  if (is.null(alpha)) {
    # type == "response": conditional mean loss rate (effective PD)
    return(stats::pnorm(stats::qnorm(p) + u_contrib))
  }

  q_alpha <- stats::qnorm(alpha)
  sq_rho <- sqrt(rho)
  sq_1mrho <- sqrt(1 - rho)

  link_mat <- outer(
    stats::qnorm(p) + u_contrib,
    sq_rho * q_alpha,
    "+"
  ) / sq_1mrho

  colnames(link_mat) <- paste0("alpha_", alpha)

  pred <- switch(type,
    link     = link_mat,
    response = stats::pnorm(link_mat)
  )

  if (ncol(pred) == 1L) {
    pred <- drop(pred)
  }
  pred
}
