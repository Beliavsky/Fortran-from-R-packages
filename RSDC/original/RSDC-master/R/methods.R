# S3 class "rsdc_fit": constructor, inference helpers, and standard methods.
# Added in 1.3-0 (P1 class/methods, P2 standard errors, P3 simulate method).

# ---- internal helpers -------------------------------------------------------

#' Observed-information variance-covariance at the MLE
#'
#' Numerically differentiates the negative log-likelihood at the final parameter
#' vector with \code{stats::optimHess} and inverts it. Returns \code{NULL} (rather
#' than erroring) when the Hessian is unavailable or singular, e.g. at a bound.
#'
#' @param par Numeric vector of final (reordered) parameters in the objective's layout.
#' @param fn Objective returning the \emph{negative} log-likelihood.
#' @param ... Passed to \code{fn}.
#' @return A variance-covariance matrix, or \code{NULL}.
#' @importFrom stats optimHess
#' @noRd
.rsdc_vcov <- function(par, fn, ndeps = 1e-3, ...) {
  # The objectives return a large finite penalty outside the feasible set, and
  # optimHess differences each coordinate by +/- ndeps. At a boundary optimum -
  # a correlation within ~2 * ndeps of +/-1, or a nearly singular correlation
  # matrix - the penalty wall falls inside the differencing step and produces a
  # huge but perfectly finite Hessian, whose inverse yields absurd (near-zero)
  # standard errors.
  #
  # The test below is the exact condition rather than a magnitude heuristic:
  # probe the objective at the points optimHess will actually visit and refuse
  # only if one of them is penalized. A magnitude threshold on the Hessian
  # cannot work here, because |H| is O(T) and O((1 - rho^2)^-2), so a
  # legitimate long-sample, high-correlation fit can exceed any fixed bound.
  f0 <- tryCatch(fn(par, ...), error = function(e) NA_real_)
  if (!is.finite(f0) || f0 >= .RSDC_PENALTY_MIN) {
    warning("Standard errors are unavailable: the optimum lies on (or outside) ",
            "the feasible boundary, where the likelihood is penalized. ",
            "Refit from a different start or interpret the fit without SEs.",
            call. = FALSE)
    return(NULL)
  }
  probes <- unlist(lapply(seq_along(par), function(i) {
    step <- c(ndeps, 2 * ndeps)
    vapply(c(par[i] + step, par[i] - step), function(v) {
      p <- par; p[i] <- v
      tryCatch(fn(p, ...), error = function(e) NA_real_)
    }, numeric(1))
  }))
  if (any(!is.finite(probes)) || any(probes >= .RSDC_PENALTY_MIN)) {
    warning("Standard errors are unavailable: the finite-difference step used ",
            "for the numerical Hessian reaches outside the feasible set (a ",
            "parameter sits within ", 2 * ndeps, " of its boundary), so the ",
            "Hessian would be penalty-contaminated. Refit at an interior optimum.",
            call. = FALSE)
    return(NULL)
  }
  H <- tryCatch(stats::optimHess(par, fn, ..., control = list(ndeps = rep(ndeps, length(par)))),
                error = function(e) NULL)
  if (is.null(H) || any(!is.finite(H))) return(NULL)
  V <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(V) || any(!is.finite(V))) return(NULL)
  V
}

#' Parameter labels matching the objective's packing order
#' @noRd
.rsdc_labels <- function(method, N, K, p) {
  pairs <- utils::combn(K, 2)
  pair_tag <- apply(pairs, 2, function(z) paste0(z[1], z[2]))
  n_reg <- if (method == "const") 1L else N
  rho_labels <- unlist(lapply(seq_len(n_reg), function(i)
    paste0("rho[g", i, ",", pair_tag, "]")))

  if (method == "const") return(rho_labels)

  if (method == "noX") {
    trans_labels <- unlist(lapply(seq_len(N), function(i) {
      if (N == 2L) paste0("p[g", i, ",stay]")
      else paste0("p[g", i, ",to", seq_len(N - 1L), "]")
    }))
    return(c(trans_labels, rho_labels))
  }

  # tvtp
  beta_labels <- unlist(lapply(seq_len(N), function(i) {
    if (N == 2L) paste0("b[g", i, ",x", seq_len(p), "]")
    else as.vector(t(outer(seq_len(N - 1L), seq_len(p),
                           function(k, x) paste0("b[g", i, ",c", k, ",x", x, "]"))))
  }))
  c(beta_labels, rho_labels)
}

#' Attach metadata, names, standard errors, and class to a backend fit
#' @noRd
.rsdc_finalize <- function(fit, method, N, K, p, call, se_requested = TRUE) {
  fit$method <- method
  fit$N      <- if (method == "const") 1L else as.integer(N)
  fit$K      <- as.integer(K)
  fit$p      <- if (method == "tvtp") as.integer(p) else NA_integer_
  fit$call   <- call

  labs <- .rsdc_labels(method, fit$N, K, p)
  if (!is.null(fit$par) && length(fit$par) == length(labs)) {
    names(fit$par) <- labs
    fit$coefficients <- fit$par
  } else {
    fit$coefficients <- fit$par
  }

  if (!is.null(fit$vcov) && all(dim(fit$vcov) == length(labs))) {
    dimnames(fit$vcov) <- list(labs, labs)
    se <- .rsdc_safe_sqrt(diag(fit$vcov))
    names(se) <- labs
    fit$se <- se
  } else {
    fit$vcov <- NULL
    fit$se   <- NULL
  }

  # Mean covariate (tvtp) used for natural-scale transition SEs and diagnostics.
  if (method == "tvtp" && !is.null(fit$X)) fit$avg_X <- colMeans(fit$X)

  # Per-observation scores for OPG / sandwich standard errors (item 8). Stored
  # whenever SEs were requested, independently of whether the Hessian inverted:
  # OPG needs only the scores, and a merely singular Hessian is the case where
  # that fallback is wanted. Excluded, however, is convergence == 99L, which
  # f_optim_noX sets when a transition row had to be projected onto the simplex:
  # that point is not a stationary optimum, so no covariance of any flavour is a
  # valid basis for inference there and the suppression must not be bypassed.
  if (isTRUE(se_requested) && !is.null(fit$residuals) &&
      !identical(as.integer(fit$convergence), 99L))
    fit$scores <- .rsdc_scores(fit$par, method, fit$residuals, fit$X, K, fit$N)

  # Filtered/smoothed regime probabilities for plot() and diagnostics (items 4, 5).
  if (method != "const" && !is.null(fit$residuals)) {
    hf <- tryCatch(rsdc_hamilton(
      y = fit$residuals,
      X = if (method == "tvtp") fit$X else NULL,
      beta = if (method == "tvtp") fit$beta else NULL,
      rho_matrix = fit$correlations, K = K, N = fit$N,
      P = if (method == "noX") fit$transition_matrix else NULL),
      error = function(e) NULL)
    if (!is.null(hf)) {
      fit$filtered_probs <- hf$filtered_probs
      fit$smoothed_probs <- hf$smoothed_probs
    }
  }

  class(fit) <- "rsdc_fit"
  fit
}

# ---- methods ---------------------------------------------------------------

#' Methods for fitted RSDC models
#'
#' A model fitted with \code{\link{rsdc_estimate}} is returned as an object of
#' class \code{"rsdc_fit"}, which supports the standard fitted-model generics —
#' so it can be used like any other fitted model in R (\code{lm}, \code{glm},
#' \dots): print and summarize it, extract coefficients, standard errors and
#' confidence intervals, compute \code{AIC}/\code{BIC}, forecast, simulate new
#' data from the fitted process, and plot the regime probabilities. Each method
#' is described in the \emph{Functions} section below.
#'
#' @param x,object An object of class \code{"rsdc_fit"}, as returned by
#'   \code{\link{rsdc_estimate}}.
#' @param digits Number of significant digits for printing.
#' @param ... Further arguments passed to methods (\code{vcov} forwards \code{B},
#'   \code{X}, \code{seed}, and \code{cores} to \code{\link{rsdc_bootstrap}} when
#'   \code{type = "bootstrap"}; otherwise unused).
#'
#' @returns \code{print} and \code{plot} return \code{x} invisibly. \code{coef},
#'   \code{nobs}, \code{logLik}, \code{vcov} and \code{confint} return the usual
#'   objects of those generics. \code{summary} returns a \code{"summary.rsdc_fit"}
#'   list with a coefficient table (estimate, SE, z, p) and regime diagnostics.
#'   \code{predict} returns the result of \code{\link{rsdc_forecast}}.
#'   \code{simulate} returns a list with \code{states} and \code{observations}.
#'
#' @examples
#' \donttest{
#' # Two persistent regimes: low (0.1) vs high (0.8) correlation
#' sim <- rsdc_simulate(n = 500, X = matrix(1, 500, 1),
#'                      beta = matrix(qlogis(0.9), 2, 1),
#'                      mu = matrix(0, 2, 2),
#'                      sigma = array(c(1, 0.1, 0.1, 1,
#'                                      1, 0.8, 0.8, 1), c(2, 2, 2)),
#'                      N = 2, seed = 2)
#' fit <- rsdc_estimate("noX", residuals = sim$observations, N = 2)
#' summary(fit)
#' confint(fit)
#' c(AIC = AIC(fit), BIC = BIC(fit))
#' plot(fit)
#' }
#'
#' @seealso \code{\link{rsdc_estimate}} for fitting;
#'   \code{\link{rsdc_forecast}}, \code{\link{rsdc_bootstrap}}.
#'
#' @name rsdc_fit-methods
NULL

#' @describeIn rsdc_fit-methods Compact printer for a fitted model.
#' @exportS3Method print rsdc_fit
print.rsdc_fit <- function(x, digits = 4, ...) {
  cat(sprintf("RSDC fit: method = \"%s\", N = %d regime(s), K = %d series, T = %d\n",
              x$method, x$N, x$K, if (is.null(x$nobs)) NA_integer_ else x$nobs))
  if (!is.null(x$call)) { cat("Call: "); print(x$call) }
  cat(sprintf("logLik = %.3f", x$log_likelihood))
  if (!is.null(x$npar)) cat(sprintf("   npar = %d   AIC = %.2f   BIC = %.2f",
                                    x$npar, stats::AIC(x), stats::BIC(x)))
  if (x$method == "const" || x$N == 1L) {
    cat("\n\nCorrelations (single regime, lower-triangular):\n")
  } else {
    cat("\n\nRegime correlations (rows = regimes, ascending mean):\n")
  }
  print(round(x$correlations, digits))
  if (x$method != "const") {
    if (x$method == "tvtp") {
      cat("\nTransition matrix (evaluated at the in-sample covariate means):\n")
    } else {
      cat("\nTransition matrix:\n")
    }
    print(round(x$transition_matrix, digits))
  }
  if (!is.null(x$convergence))
    cat(sprintf("\nOptimiser convergence code: %d (0 = converged)\n", x$convergence))
  if (is.null(x$se))
    cat("Standard errors: not available (compute_se = FALSE or singular/unavailable Hessian).\n")
  invisible(x)
}

#' @describeIn rsdc_fit-methods Coefficient vector (named, in objective order).
#' @exportS3Method coef rsdc_fit
coef.rsdc_fit <- function(object, ...) object$coefficients

#' @describeIn rsdc_fit-methods Number of observations used in estimation.
#' @exportS3Method nobs rsdc_fit
nobs.rsdc_fit <- function(object, ...) object$nobs

#' @describeIn rsdc_fit-methods Log-likelihood (carries \code{df} and \code{nobs} so
#'   \code{stats::AIC}/\code{BIC} work out of the box).
#' @exportS3Method logLik rsdc_fit
logLik.rsdc_fit <- function(object, ...) {
  val <- object$log_likelihood
  attr(val, "df")    <- object$npar
  attr(val, "nobs")  <- object$nobs
  class(val) <- "logLik"
  val
}

#' @describeIn rsdc_fit-methods Variance-covariance matrix of the estimates. The
#'   \emph{numerical} estimators are \code{type = "hessian"} (default, inverse
#'   observed information), \code{"opg"} (outer product of gradients) and
#'   \code{"sandwich"} (QML/robust \eqn{H^{-1} (\sum_t s_t s_t') H^{-1}}); the
#'   \emph{simulation-based} estimator is \code{"bootstrap"}, which calls
#'   \code{\link{rsdc_bootstrap}} (pass \code{B}, \code{X}, \code{seed}, and
#'   \code{cores} through \code{...}). The bootstrap is recomputed on each call.
#' @param type Covariance estimator: one of \code{"hessian"}, \code{"opg"},
#'   \code{"sandwich"} (numerical) or \code{"bootstrap"} (simulation-based).
#' @exportS3Method vcov rsdc_fit
vcov.rsdc_fit <- function(object, type = c("hessian", "opg", "sandwich", "bootstrap"), ...) {
  type <- match.arg(type)
  if (type == "bootstrap") {
    dots  <- list(...)
    B     <- if (!is.null(dots$B)) dots$B else 199L
    cores <- if (!is.null(dots$cores)) dots$cores else 1
    return(rsdc_bootstrap(object, B = B, X = dots$X, seed = dots$seed,
                          cores = cores)$vcov)
  }
  # "hessian" and "sandwich" need the Hessian bread; "opg" needs only the
  # per-observation scores, so it remains available when the Hessian alone is
  # singular (scores are withheld where no covariance is valid, e.g. a projected
  # non-stationary point, and .rsdc_robust_vcov then errors informatively).
  if (is.null(object$vcov) && !(type == "opg" && !is.null(object$scores)))
    stop("No variance-covariance available (Hessian singular/unavailable, or ",
         "the fit did not converge to an interior optimum). Refit with ",
         "control = list(compute_se = TRUE) from a different start.")
  if (is.null(object$vcov) && type == "opg" &&
      !identical(as.integer(object$convergence), 0L))
    warning("Using OPG standard errors at a fit that did not converge cleanly ",
            "(convergence = ", object$convergence, "); interpret with caution.",
            call. = FALSE)
  V <- .rsdc_robust_vcov(object, type)
  if (is.null(V)) stop("Requested covariance (", type, ") is unavailable (singular).")
  V
}

#' @describeIn rsdc_fit-methods Wald confidence intervals from the chosen covariance
#'   (numerical or \code{"bootstrap"}). For bootstrap \emph{percentile} intervals
#'   instead of Wald, use \code{\link{rsdc_bootstrap}} directly.
#' @param parm Vector of parameter names/indices (default: all).
#' @param level Confidence level.
#' @exportS3Method confint rsdc_fit
#' @importFrom stats qnorm
confint.rsdc_fit <- function(object, parm, level = 0.95,
                             type = c("hessian", "opg", "sandwich", "bootstrap"), ...) {
  type <- match.arg(type)
  cf <- coef(object); V <- vcov(object, type = type, ...)
  ses <- .rsdc_safe_sqrt(diag(V))  # consistent with summary(); non-PD diagonals -> NA
  if (missing(parm)) parm <- names(cf)
  a <- (1 - level) / 2
  fac <- stats::qnorm(c(a, 1 - a))
  ci <- cf[parm] + outer(ses[parm], fac)
  colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE), "%")
  rownames(ci) <- parm
  ci
}

#' @describeIn rsdc_fit-methods Summary with a coefficient table (estimate, SE, z, p).
#' @exportS3Method summary rsdc_fit
#' @importFrom stats pnorm
summary.rsdc_fit <- function(object, ...) {
  cf <- coef(object)
  tab <- NULL
  if (!is.null(object$se)) {
    se <- object$se[names(cf)]
    z  <- cf / se
    pv <- 2 * stats::pnorm(-abs(z))
    tab <- cbind(Estimate = cf, `Std. Error` = se, `z value` = z, `Pr(>|z|)` = pv)
  }
  out <- list(call = object$call, method = object$method, N = object$N, K = object$K,
              nobs = object$nobs, npar = object$npar,
              logLik = object$log_likelihood,
              AIC = stats::AIC(object), BIC = stats::BIC(object),
              coefficients = tab, correlations = object$correlations,
              transition_matrix = object$transition_matrix, convergence = object$convergence,
              natural_se = .rsdc_natural_se(object),
              diagnostics = if (object$method != "const") .rsdc_diagnostics(object) else NULL)
  class(out) <- "summary.rsdc_fit"
  out
}

#' @exportS3Method print summary.rsdc_fit
#' @importFrom stats printCoefmat
print.summary.rsdc_fit <- function(x, digits = 4, ...) {
  cat(sprintf("RSDC fit summary: method = \"%s\", N = %d, K = %d, T = %s\n",
              x$method, x$N, x$K, ifelse(is.null(x$nobs), "NA", x$nobs)))
  cat(sprintf("logLik = %.3f   npar = %s   AIC = %.2f   BIC = %.2f\n\n",
              x$logLik, ifelse(is.null(x$npar), "NA", x$npar), x$AIC, x$BIC))
  if (!is.null(x$coefficients)) {
    cat("Coefficients:\n")
    stats::printCoefmat(x$coefficients, digits = digits, has.Pvalue = TRUE)
  } else {
    cat("Coefficients: standard errors unavailable (compute_se = FALSE or singular Hessian).\n")
    cat("Regime correlations:\n")
    print(round(x$correlations, digits))
  }
  if (!is.null(x$natural_se) && !is.null(x$natural_se$transition)) {
    cat(if (identical(x$method, "tvtp")) {
      "\nTransition probabilities at the in-sample covariate means (delta-method SE):\n"
    } else {
      "\nTransition probabilities (delta-method SE):\n"
    })
    print(format(x$natural_se$transition, digits = digits))
  }
  if (!is.null(x$diagnostics)) {
    cat(if (identical(x$method, "tvtp")) {
      "\nRegime diagnostics at the in-sample covariate means (stay prob., expected duration, ergodic prob.):\n"
    } else {
      "\nRegime diagnostics (stay prob., expected duration, ergodic prob.):\n"
    })
    print(format(x$diagnostics, digits = digits))
  }
  invisible(x)
}

#' @describeIn rsdc_fit-methods Forecast from a fitted model (wraps
#'   \code{\link{rsdc_forecast}}); supply the residuals, conditional volatilities,
#'   and (for \code{"tvtp"}) the covariate matrix \code{X}.
#' @param residuals Numeric matrix \eqn{T \times K} of standardized residuals
#'   on which the forecast filter is run.
#' @param sigma_matrix Numeric matrix of conditional standard deviations.
#' @param value_cols Columns of \code{sigma_matrix} giving the asset order.
#' @param X Covariate matrix \eqn{T \times p} (required for a \code{"tvtp"} fit;
#'   \code{simulate} uses its row count as the simulation length).
#' @param out_of_sample Logical; passed to \code{\link{rsdc_forecast}}
#'   (70/30 split).
#' @exportS3Method predict rsdc_fit
predict.rsdc_fit <- function(object, residuals, sigma_matrix, value_cols,
                             X = NULL, out_of_sample = FALSE, ...) {
  rsdc_forecast(method = object$method, N = object$N, residuals = residuals,
                X = X, final_params = object, sigma_matrix = sigma_matrix,
                value_cols = value_cols, out_of_sample = out_of_sample, ...)
}

#' @describeIn rsdc_fit-methods Simulate from a fitted model. For \code{"tvtp"} supply
#'   a covariate matrix \code{X} (its row count sets the length); for \code{"noX"}
#'   the fixed transition matrix is used; for \code{"const"} a single regime is drawn.
#' @param nsim Unused (kept for generic compatibility; one path is returned).
#' @param seed Optional RNG seed.
#' @param n Series length for \code{"noX"}/\code{"const"} (defaults to the fitted T).
#' @exportS3Method simulate rsdc_fit
simulate.rsdc_fit <- function(object, nsim = 1, seed = NULL, X = NULL, n = NULL, ...) {
  if (!is.null(nsim) && !isTRUE(nsim == 1))
    warning("simulate(): nsim > 1 is not supported; returning a single simulated path.")
  if (!is.null(seed)) set.seed(seed)
  N <- object$N; K <- object$K
  sigma <- object$covariances
  mu <- matrix(0, max(N, 1L), K)

  if (object$method == "tvtp") {
    if (is.null(X)) stop("simulate(): provide X (n x p) for a 'tvtp' fit.")
    return(rsdc_simulate(n = nrow(X), X = X, beta = object$beta,
                         mu = mu, sigma = sigma, N = N, seed = NULL))
  }

  n_ <- if (!is.null(n)) as.integer(n) else if (!is.null(object$nobs)) object$nobs else 500L
  P <- object$transition_matrix
  states <- integer(n_); obs <- matrix(NA_real_, n_, K)
  # Draw the initial noX state from the chain's ergodic distribution (uniform fallback).
  if (object$method == "noX") {
    pi0 <- .rsdc_ergodic(P)
    if (any(!is.finite(pi0))) pi0 <- rep(1 / N, N)
    states[1] <- sample.int(N, 1, prob = pi0)
  } else {
    states[1] <- 1L
  }
  obs[1, ] <- mvtnorm::rmvnorm(1, mu[states[1], ], sigma[, , states[1]])
  if (n_ > 1L) for (t in 2:n_) {
    states[t] <- if (object$method == "noX") sample.int(N, 1, prob = P[states[t - 1L], ]) else 1L
    obs[t, ] <- mvtnorm::rmvnorm(1, mu[states[t], ], sigma[, , states[t]])
  }
  list(states = states, observations = obs, transition_matrices = NULL)
}

#' @describeIn rsdc_fit-methods Plot the smoothed regime probabilities (one panel per
#'   regime) for a fitted \code{"noX"}/\code{"tvtp"} model. \code{which = "filtered"}
#'   plots filtered probabilities instead.
#' @param which Either \code{"smoothed"} (default) or \code{"filtered"}.
#' @exportS3Method plot rsdc_fit
#' @importFrom graphics par plot abline
plot.rsdc_fit <- function(x, which = c("smoothed", "filtered"), ...) {
  which <- match.arg(which)
  pr <- if (which == "filtered") x$filtered_probs else x$smoothed_probs
  if (is.null(pr))
    stop("No regime probabilities stored. The constant model has a single regime; ",
         "otherwise the fit lacks stored residuals or the filter failed at the ",
         "estimates. Use rsdc_hamilton() on your residuals to compute them directly.")
  N <- nrow(pr); Tn <- ncol(pr)
  op <- graphics::par(mfrow = c(N, 1), mar = c(3.2, 4, 1.4, 1), las = 1)
  on.exit(graphics::par(op))
  for (i in seq_len(N)) {
    graphics::plot(seq_len(Tn), pr[i, ], type = "l", ylim = c(0, 1),
                   xlab = if (i == N) "Time" else "", ylab = sprintf("P(regime %d)", i),
                   main = sprintf("%s regime-%d probability", which, i), ...)
    graphics::abline(h = c(0, 1), col = "grey80", lty = 3)
  }
  invisible(x)
}
