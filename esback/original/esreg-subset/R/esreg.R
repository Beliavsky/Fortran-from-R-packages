#' Joint Quantile and Expected Shortfall Regression
#'
#' Estimates a joint linear regression model for the pair (VaR, ES):
#' \deqn{Q_\alpha(Y | Xq) = Xq'\beta_q}
#' \deqn{ES_\alpha(Y | Xe) = Xe'\beta_e}
#'
#' @param formula Formula: y ~ x1 + x2 ... | x1 + x2 ...
#' where the first part after the response variable specifies the quantile equation
#' and the second the expected shortfall part. If only one set of regressors is
#' provided it is used for both model specifications.
#' @param data data.frame that holds the variables
#' @param y Response vector
#' @param xq Explanatory variables for the quantile regression equation
#' @param xe Explanatory variables for the expected shortfall regression equation
#' @param alpha Probability level
#' @param g1 1, 2 (see \link{G1_fun}, \link{G1_prime_fun}), defaults to 1
#' @param g2 1, 2, 3, 4, 5 (see \link{G2_curly_fun}, \link{G2_fun}, \link{G2_prime_fun}). defaults to 2
#' @param early_stopping Stop the iterated local search if there is no improvement in early_stopping steps.
#' @param ... Further arguments (does not apply here)
#' @return An esreg object
#' @seealso \code{\link{vcov.esreg}} for covariance estimation
#' @examples
#' # Simulate data (DGP-(2) in the linked paper)
#' set.seed(0)
#' x <- rchisq(1000, df=1)
#' y <- -x + (1 + 0.5 * x) * rnorm(1000)
#'
#' # True quantile and expected shortfall regression parameters (for alpha=0.025)
#' alpha=0.025
#' true_pars <- c(-1.959964, -1.979982, -2.337803, -2.168901)
#'
#' # Estimate the model using the standard settings
#' fit <- esreg(y ~ x, alpha=alpha)
#'
#' # Compare the different variance-covariance estimators
#' cov1 <- vcov(object=fit, sparsity="iid", sigma_est="ind")
#' cov2 <- vcov(object=fit, sparsity="nid", sigma_est="scl_N")
#' cov3 <- vcov(object=fit, sparsity="nid", sigma_est="scl_sp")
#'
#' print("Comparison of the variance-covariance estimators")
#' print(cbind(Truth=true_pars,
#'             Estimate=coef(fit),
#'             SE_iid_ind=sqrt(diag(cov1)),
#'             SE_nid_N=sqrt(diag(cov2)),
#'             SE_nid_sp=sqrt(diag(cov3))))
#'
#' # Compares estimates using different G2 functions
#' fit1 <- esreg(y ~ x, alpha=alpha, g2=1)
#' fit2 <- esreg(y ~ x, alpha=alpha, g2=2)
#' fit3 <- esreg(y ~ x, alpha=alpha, g2=3)
#' fit4 <- esreg(y ~ x, alpha=alpha, g2=4)
#' fit5 <- esreg(y ~ x, alpha=alpha, g2=5)
#' fits <- sapply(list(fit1, fit2, fit3, fit4, fit5), coef)
#' colnames(fits) <- sapply(1:5, function(i) esreg:::.G_function_names(1, i)[2])
#' print("Comparison of the five G2 functions")
#' print(rbind(Truth=true_pars, t(fits)))
#'
#' # Usage of different covariates
#' x <- rchisq(1000, df=1)
#' noise <- rnorm(1000)
#' y <- -x + (1 + 0.5 * x) * rnorm(1000)
#' fit <- esreg(y ~ x | x + noise, alpha=0.025)
#' print("Using different covariates for VaR and ES")
#' print(summary(fit))
#'
#' @references \href{https://arxiv.org/abs/1704.02213}{A Joint Quantile and Expected Shortfall Regression Framework}
#' @rdname esreg
#' @export
esreg <- function(...) UseMethod('esreg')

#' @rdname esreg
#' @export
esreg.formula <- function(formula, data=parent.frame(), alpha, g1 = 2L, g2 = 1L,
                          early_stopping = 10, ...) {
  chkDots(...)

  # Prepare the formula
  if(is.matrix(data)) {
    data <- as.data.frame(data)
  }
  formula <- Formula::Formula(formula)
  if (length(formula)[2] == 1) {
    formula <- Formula::as.Formula(formula(formula),
                                   formula(formula, lhs = 0, rhs = 1))
  } else if (length(formula)[2] > 2) {
    stop('You must not have more than two formula objects.')
  }
  mf <- stats::model.frame(formula = formula, data = data, na.action = stats::na.pass)
  terms <- stats::terms(formula)

  # Extract the data matrices
  xq <- stats::model.matrix(formula, mf, rhs=1)
  xe <- stats::model.matrix(formula, mf, rhs=2)
  y <- stats::model.response(mf)

  # Fit the model
  fit <- esreg.fit(xq, xe, y, alpha, g1, g2, early_stopping)
  fit$call <- match.call()
  fit$terms <- terms
  fit$formula <- formula
  fit$data <- data

  fit
}

#' @rdname esreg
#' @export
esreg.default <- function(xq, xe, y, alpha, g1 = 2L, g2 = 1L,
                          early_stopping = 10, ...) {
  chkDots(...)

  xq <- cbind("(Intercept)" = 1, xq)
  xe <- cbind("(Intercept)" = 1, xe)

  # Fit the model
  fit <- esreg.fit(xq, xe, y, alpha, g1, g2, early_stopping)
  fit$call <- match.call()

  fit
}

#' @export
model.matrix.esreg <- function(object, ...) {
  chkDots(...)
  mm1 <- stats::model.matrix(object$formula, rhs=1, data=object$data)
  mm2 <- stats::model.matrix(object$formula, rhs=2, data=object$data)
  mm <- cbind(mm1, mm2)
  attributes(mm)$assign <- c(attr(mm1,  'assign'), attr(mm2,  'assign'))

  mm
}

#' @export
print.esreg <- function(x, digits = 4, ...) {
  chkDots(...)
  cat("Call:\n")
  cat(deparse(x$call), "\n")
  cat("\nQuantile Coefficients:\n")
  print(format(x$coefficients_q, digits = digits), quote = FALSE)
  cat("\nExpected Shortfall Coefficients:\n")
  print(format(x$coefficients_e, digits = digits), quote = FALSE)
}

#' @export
summary.esreg <- function(object, ...) {
  chkDots(...)
  xq <- object$xq
  xe <- object$xe
  n <- nrow(xq)
  k <- ncol(xq) + ncol(xe)
  cov <- vcov.esreg(object, ...)
  se <- sqrt(diag(cov))
  tval <- stats::coef(object) / se
  coef_mat <- cbind(
    Estimate     = stats::coef(object),
    `Std. Error` = se,
    `t value`    = tval,
    `Pr(>|t|)`   = 2 * stats::pt(abs(tval), n - 2*k, lower.tail = FALSE)
  )

  structure(
    c(object, list(cov = cov, coef_mat = coef_mat)),
    class = "summary.esreg")
}

#' @export
print.summary.esreg <- function(x, ...) {
  kq <- length(x$coefficients_q)
  ke <- length(x$coefficients_e)
  cat("Call:\n", paste0(deparse(x$call), sep = "\n", collapse = "\n"))
  cat("\nQuantile Coefficients:\n")
  stats::printCoefmat(x$coef_mat[1:kq,,drop=FALSE], signif.legend = FALSE)
  cat("\nExpected Shortfall Coefficients:\n")
  stats::printCoefmat(x$coef_mat[(kq+1):(kq+ke),,drop=FALSE])
}


#' @export
fitted.esreg <- function(object, ...) {
  chkDots(...)
  cbind(object$xq %*% object$coefficients_q,
        object$xe %*% object$coefficients_e)
}

#' @export
residuals.esreg <- function(object, ...) {
  chkDots(...)
  cbind(object$y - object$xq %*% object$coefficients_q,
        object$y - object$xe %*% object$coefficients_e)
}

#' @export
predict.esreg <- function(object, newdata=NULL, ...) {
  chkDots(...)
  if (is.null(newdata)) {
    yhat <- fitted.esreg(object)
  } else {
    if (nrow(newdata) == 0)
      stop('newdata is empty')
    if (is.null(object$formula)) {
      stop('The predict method is only supported using the formula interface')
    } else {
      mf <- stats::model.frame(object$formula, newdata)
      xq <- stats::model.matrix(object$formula, mf, rhs=1)
      xe <- stats::model.matrix(object$formula, mf, rhs=2)
      yhat <- cbind(xq %*% object$coefficients_q,
                    xe %*% object$coefficients_e)
    }
  }

  yhat
}


#' Estimating function
#'
#' This function matches the estfun function of the sandwich package and
#' returns the estimating functions for the fitted model.
#' It can for instance be used for an OPG estimator of the sigma matrix.
#' For esreg, the dimension of the estimating functions is n x (kq + ke).
#'
#' @param x An \link{esreg} object
#' @param ... Further arguments (does not apply here)
#' @export
estfun.esreg <- function(x, ...) {
  chkDots(...)

  # Extract elements from object
  object <- x
  y <- object$y
  xq <- object$xq
  xe <- object$xe
  coefficients_q <- object$coefficients_q
  coefficients_e <- object$coefficients_e
  alpha <- object$alpha
  n <- length(y)
  kq <- ncol(xq)
  ke <- ncol(xe)

  # Transform the data and coefficients
  if (object$g2 %in% c(1, 2, 3)) {
    max_y <- max(y)
    y <- y - max_y
    coefficients_q[1] <- coefficients_q[1] - max_y
    coefficients_e[1] <- coefficients_e[1] - max_y
  }

  # Precompute some quantities
  xbq <- as.numeric(xq %*% coefficients_q)
  xbe <- as.numeric(xe %*% coefficients_e)
  uq <- as.numeric(y - xbq)

  # Evaluate G1 / G2 functions
  G1_prime_xq <- G_vec(z = xbq, g = "G1_prime", type = object$g1)
  G2_xe <- G_vec(z = xbe, g = "G2", type = object$g2)
  G2_prime_xe <- G_vec(z = xbe, g = "G2_prime", type = object$g2)
  G2_prime_prime_xe <- G_vec(z = xbe, g = "G2_prime_prime", type = object$g2)

  psi <- estimating_function_loop(
    y = y, xq = xq, xe = xe, xbq = xbq, xbe = xbe, alpha = alpha,
    G1_prime_xq = G1_prime_xq,
    G2_xe = G2_xe, G2_prime_xe = G2_prime_xe)

  psi
}

#' Covariance Estimation
#'
#' Estimate the variance-covariance matrix of the joint (VaR, ES) estimator
#'
#' @param object An \link{esreg} object
#' @param method For asymptotic use \link{vcovA}, for boot use \link{vcovB}
#' @param ... All possible values which can be passed to \link{vcovA} and \link{vcovB}
#' @export
vcov.esreg <- function(object, method='asymptotic', ...) {
  if (method == 'asymptotic') {
    cov <- vcovA(object, ...)
  } else if(method == 'boot') {
    cov <- vcovB(object, ...)
  } else {
    stop('method can be asymptotic or boot')
  }
  cov
}

#' Asymptotic Covariance Estimation
#'
#' Estimate the variance-covariance matrix of the joint (VaR, ES) estimator by the sandwich formula:
#' \deqn{\lambda^{-1} \Sigma \lambda^{-1}}
#' Several estimators are available for both matrices and the default options are selected to take into account
#' possible misspecifications in the underlying data.
#'
#' @param object An esreg object
#' @param sigma_est The estimator to be used for \eqn{\Sigma}, see \link{conditional_truncated_variance}
#'   \itemize{
#'     \item ind - Variance over all negative residuals
#'     \item scl_N - Scaling with the normal distribution
#'     \item scl_sp - Scaling with the kernel density function
#'     }
#' @param sparsity The estimator to be used for the sparsity in \eqn{\Lambda}, see \link{density_quantile_function}
#'   \itemize{
#'     \item iid - Piecewise linear interpolation of the distribution
#'     \item nid - Hendricks and Koenker sandwich
#'   }
#' @param bandwidth_estimator The bandwidth estimator to be used for the iid and nid sparsity estimator, see \link{density_quantile_function}
#'  \itemize{
#'    \item Bofinger
#'    \item Chamberlain
#'    \item Hall-Sheather
#'  }
#' @param misspec if TRUE, the estimator accounts for potential misspecification in the model
#' @export
vcovA <- function(object, sigma_est = 'scl_sp', sparsity = 'nid', misspec = TRUE, bandwidth_estimator = 'Hall-Sheather') {
  lambda <- lambda_matrix(object = object, sparsity = sparsity,
                          bandwidth_estimator = bandwidth_estimator, misspec = misspec)
  lambda_inverse <- solve(lambda)
  sigma <- sigma_matrix(object = object, sigma_est = sigma_est, misspec = misspec)
  n <- length(object$y)
  cov <- 1/n * (lambda_inverse %*% sigma %*% lambda_inverse)
  rownames(cov) <- colnames(cov) <- names(stats::coef(object))
  cov
}

#' @title Bootstrap Covariance Estimation
#' @description Estimate the variance-covariance matrix of the joint (VaR, ES) estimator
#' using the bootstrap.
#' @param object An esreg object
#' @param bootstrap_method The bootstrap sampling scheme to be used
#'   \itemize{
#'     \item iid - The iid bootstrap of Efron (1979)
#'   }
#' @param B The number of bootstrap iterations
#' @export
vcovB <-function(object, bootstrap_method='iid', B=1000) {
  if (!(bootstrap_method %in% c("iid")))
    stop("bootstrap_method can be iid")
  if (B < 1000)
    warning("The number of bootstrap iterations is small!")

  # Draw the bootstrap indices
  n <- length(object$y)
  idx <- matrix(sample(1:n, size = n * B, replace = TRUE), nrow = n)

  # Get the bootstrap parameter estimates
  b <- sapply(seq_len(B), function(i) {
    fitb <- esreg(object$y[idx[,i]] ~ object$xq[idx[,i], -1] | object$xe[idx[,i], -1],
                  alpha = object$alpha, g1 = object$g1, g2 = object$g2,
                  early_stopping = 0)
    fitb$coefficients
  })

  # Compute the covariance
  cov <- cov(t(b))
  rownames(cov) <- colnames(cov) <- names(stats::coef(object))
  cov
}

#' @export
esreg.fit <- function(xq, xe, y, alpha, g1, g2, early_stopping, ...) {
  t0 <- Sys.time()

  # Check input parameters and data
  if (!(g1 %in% c(1, 2))) stop("G1 can be 1 or 2.")
  if (!(g2 %in% c(1, 2, 3, 4, 5))) stop("G2 can be 1, 2, 3, 4 or 5.")
  if ((alpha < 0) | (1 < alpha)) stop("alpha not in (0, 1)")
  if (any(is.na(y)) | any(is.na(xq)) | any(is.na(xe)))
    stop("Data contains NAs!")
  if (!(all(is.finite(y)) & all(is.finite(xq)) & all(is.finite(xe))))
    stop("Not all values are finite!")
  if (length(y) != nrow(xq))
    stop("Dimensions of y and xq do not match!")
  if (length(y) != nrow(xe))
    stop("Dimensions of y and xe do not match!")

  # Transform the data
  if (g2 %in% c(1, 2, 3)) {
    max_y <- max(y)
    y <- y - max_y
  }

  ## Find starting values
  # Match quantile and expected shortfall levels
  e <- -stats::dnorm(stats::qnorm(alpha))/alpha
  alpha_tilde <- stats::uniroot(function(x) stats::qnorm(x) - e, c(0, alpha))$root

  # Fit the quantile regression
  fit_qr_q <- suppressWarnings(quantreg::rq(y ~ xq - 1, tau = alpha))
  fit_qr_e <- suppressWarnings(quantreg::rq(y ~ xe - 1, tau = alpha_tilde))

  # Get standard errors for the perturbation
  qr_sum_q <- suppressWarnings(summary(fit_qr_q, se="iid"))
  qr_sum_e <- suppressWarnings(summary(fit_qr_e, se="iid"))
  se <- c(qr_sum_q$coefficients[,2], qr_sum_e$coefficients[,2])

  # Store the coefficient estimates
  b0 <- c(fit_qr_q$coef, fit_qr_e$coef)

  # If G2 is 1, 2, 3 ensure that x'be < 0 by moving the es intercept
  kq <- ncol(xq)
  ke <- ncol(xe)
  if (g2 %in% c(1, 2, 3)) {
    max_xe <- max(xe %*% b0[(kq+1):(kq+ke)])
    b0[kq+1] <- b0[kq+1] - (max_xe + 0.1) * (max_xe >= -0.1)
  }

  # Set the target function
  fun <- function(b) suppressWarnings(
    esr_rho_lp(b = b, y = y, xq = xq, xe = xe, alpha = alpha, g1 = g1, g2 = g2)
  )

  ## Optimize the model using iterated local search
  # Initial optimization
  fit <- try(stats::optim(par = b0, fn = fun, method = "Nelder-Mead",
                          control = list(maxit = 2000)), silent = TRUE)

  # Counts the iterations without decrease of the loss
  counter <- 0
  while (counter < early_stopping) {
    # Perturbe b
    bt <- fit$par + stats::rnorm(kq+ke, sd=se)

    # If G2 is 1, 2, 3 ensure that x'be < 0 by moving the es intercept
    if (g2 %in% c(1, 2, 3)) {
      max_xe <- max(xe %*% bt[(kq+1):(kq+ke)])
      bt[kq+1] <- bt[kq+1] - (max_xe + 0.1) * (max_xe >= 0)
    }

    # Fit the model with the perturbed parameters
    tmp_fit <- try(stats::optim(par = bt, fn = fun, method = "Nelder-Mead",
                                control = list(maxit = 2000)), silent = TRUE)

    # Replace the fit if the new loss is smaller than the old. Otherwise increase the counter.
    if (!inherits(tmp_fit, "try-error")) {
      if (tmp_fit$value < fit$value) {
        fit <- tmp_fit
        counter <- 0
      } else {
        counter <- counter + 1
      }
    }
  }

  # Set names of the parameters
  names(fit$par) <- c(paste0("bq_", 1:kq - 1), paste0("be_", 1:ke - 1))

  # Undo the transformation
  if (g2 %in% c(1, 2, 3)) {
    y <- y + max_y
    fit$par[1] <- fit$par[1] + max_y
    fit$par[kq+1] <- fit$par[kq+1] + max_y
  }

  # Retun results
  revtal <- list(
    coefficients   = fit$par,
    coefficients_q = fit$par[1:kq],
    coefficients_e = fit$par[(kq+1):(kq+ke)],
    y              = y,
    xq             = xq,
    xe             = xe,
    alpha          = alpha,
    g1             = g1,
    g2             = g2,
    loss           = fit$value,
    time           = Sys.time() - t0
  )

  class(revtal) <- 'esreg'

  revtal
}

#' Lambda Matrix
#'
#' Estimate the lambda matrix.
#'
#' @inheritParams vcovA
#' @export
lambda_matrix <- function(object, sparsity, bandwidth_estimator, misspec) {
  if(!(sparsity %in% c("iid", "nid")))
    stop("sparsity can be iid or nid")
  if(!(bandwidth_estimator %in% c("Bofinger", "Chamberlain", "Hall-Sheather")))
    stop("bandwidth_estimator can be Bofinger, Chamberlain or Hall-Sheather")

  # Extract elements from object
  y <- object$y
  xq <- object$xq
  xe <- object$xe
  coefficients_q <- object$coefficients_q
  coefficients_e <- object$coefficients_e
  alpha <- object$alpha
  n <- length(y)
  kq <- ncol(xq)
  ke <- ncol(xe)

  # Transform the data and coefficients
  if (object$g2 %in% c(1, 2, 3)) {
    max_y <- max(y)
    y <- y - max_y
    coefficients_q[1] <- coefficients_q[1] - max_y
    coefficients_e[1] <- coefficients_e[1] - max_y
  }

  # Precompute some quantities
  xbq <- as.numeric(xq %*% coefficients_q)
  xbe <- as.numeric(xe %*% coefficients_e)
  uq <- as.numeric(y - xbq)

  # Evaluate G1 / G2 functions
  G1_prime_xq <- G_vec(z = xbq, g = "G1_prime", type = object$g1)
  G1_prime_prime_xq <- G_vec(z = xbq, g = "G1_prime_prime", type = object$g1)
  G2_xe <- G_vec(z = xbe, g = "G2", type = object$g2)
  G2_prime_xe <- G_vec(z = xbe, g = "G2_prime", type = object$g2)
  G2_prime_prime_xe <- G_vec(z = xbe, g = "G2_prime_prime", type = object$g2)

  # Check the methods in case of sample quantile / es
  if ((kq == 1) & (ke == 1) & sparsity != "iid") {
    warning("Changed sparsity estimation to iid!")
    sparsity <- "iid"
  }

  # Density quantile function
  dens <- density_quantile_function(y = y, x = xq, u = uq, alpha = object$alpha,
                                    sparsity = sparsity, bandwidth_estimator = bandwidth_estimator)

  # Conditional CDF evaluated at conditional quantile
  cdf <- cdf_at_quantile(y = y, x = xq, q = xbq)

  # Compute lambda
  lambda <- lambda_matrix_loop(
    xq = xq, xe = xe, xbq = xbq, xbe = xbe, alpha = alpha,
    G1_prime_xq = G1_prime_xq, G1_prime_prime_xq = G1_prime_prime_xq,
    G2_xe = G2_xe, G2_prime_xe = G2_prime_xe, G2_prime_prime_xe = G2_prime_prime_xe,
    density = dens, cdf = cdf, include_misspecification_terms = misspec)

  lambda
}

#' Sigma Matrix
#'
#' Estimate the sigma matrix.
#'
#' @inheritParams vcovA
#' @export
sigma_matrix <- function(object, sigma_est, misspec) {
  if(!(sigma_est %in% c("ind", "scl_N", "scl_sp")))
    stop("sigma_estimator can be ind, scl_N or scl_sp")

  # Extract elements from object
  y <- object$y
  xq <- object$xq
  xe <- object$xe
  coefficients_q <- object$coefficients_q
  coefficients_e <- object$coefficients_e
  alpha <- object$alpha
  n <- length(y)
  kq <- ncol(xq)
  ke <- ncol(xe)

  # Transform the data and coefficients
  if (object$g2 %in% c(1, 2, 3)) {
    max_y <- max(y)
    y <- y - max_y
    coefficients_q[1] <- coefficients_q[1] - max_y
    coefficients_e[1] <- coefficients_e[1] - max_y
  }

  # Precompute some quantities
  xbq <- as.numeric(xq %*% coefficients_q)
  xbe <- as.numeric(xe %*% coefficients_e)
  uq <- as.numeric(y - xbq)

  # Evaluate G1 / G2 functions
  G1_prime_xq <- G_vec(z = xbq, g = "G1_prime", type = object$g1)
  G2_xe <- G_vec(z = xbe, g = "G2", type = object$g2)
  G2_prime_xe <- G_vec(z = xbe, g = "G2_prime", type = object$g2)
  G2_prime_prime_xe <- G_vec(z = xbe, g = "G2_prime_prime", type = object$g2)

  # Check the methods in case of sample quantile / es
  if ((kq == 1) & (ke == 1) & sigma_est != "ind") {
    warning("Changed conditional truncated variance estimation to ind!")
    sigma_est <- "ind"
  }

  # Estimate the (conditional) truncated variance
  cv <- conditional_truncated_variance(y = uq, x = xq, approach = sigma_est)

  # Estimate the CDF at the quantile predictions
  cdf <- cdf_at_quantile(y = y, x = xq, q = xbq)

  # Compute sigma
  sigma <- sigma_matrix_loop(
    xq = xq, xe = xe, xbq = xbq, xbe = xbe, alpha = alpha,
    G1_prime_xq = G1_prime_xq,
    G2_xe = G2_xe, G2_prime_xe = G2_prime_xe,
    conditional_variance = cv, cdf = cdf,
    include_misspecification_terms = misspec)

  sigma
}

