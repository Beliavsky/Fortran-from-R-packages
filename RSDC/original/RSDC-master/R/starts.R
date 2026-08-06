#' Data-driven warm starts for high-dimensional estimation
#'
#' Builds a set of feasible starting vectors for \code{\link{rsdc_estimate}}
#' from the empirical correlations of low/…/high average-correlation
#' sub-periods. Pass the result via \code{control = list(start = )} to run a
#' warm-started multi-start search (see \emph{Details}).
#'
#' In moderate dimensions (\eqn{K \le 4}) the default global search
#' (\pkg{DEoptim}) works well. In higher dimensions it breaks down: a vector of
#' \eqn{K(K-1)/2} pairwise correlations drawn uniformly in \eqn{(-1,1)} is
#' almost never a valid (positive-definite) correlation matrix, so the random
#' initial population contains no feasible point. \code{rsdc_starts} sidesteps
#' this by initializing each regime at an \emph{empirical} correlation matrix,
#' which is positive semi-definite by construction and positive definite
#' whenever the group has enough non-collinear observations (starts that are not
#' usable are filtered out by the feasibility check below):
#' \enumerate{
#'   \item compute the rolling mean pairwise correlation over \code{window}
#'         observations;
#'   \item split the sample into \code{N} groups by quantiles of that series
#'         (calm, …, turbulent), so regimes are ordered by ascending
#'         correlation, matching the package's identification convention;
#'   \item use each group's empirical correlation matrix as that regime's
#'         starting values, with persistent initial transitions (stay
#'         probability \code{stay}) and, for \code{"tvtp"}, zero covariate
#'         slopes (the data reveal the covariate effect during estimation);
#'   \item replicate the start under \code{n_starts} shrinkage factors of the
#'         correlations towards zero (which preserves positive definiteness),
#'         spanning "regimes as contrasted as the split suggests" to "regimes
#'         much closer together" — protection against local optima.
#' }
#'
#' @param residuals Numeric matrix \eqn{T \times K} of standardized residuals,
#'   as passed to \code{\link{rsdc_estimate}}.
#' @param N Integer. Number of regimes (ignored for \code{method = "const"}).
#' @param method Character. One of \code{"tvtp"}, \code{"noX"}, \code{"const"}.
#' @param X Covariate matrix \eqn{T \times p} (required for \code{"tvtp"}).
#' @param window Integer. Rolling-correlation window used to classify
#'   observations into regimes (default 126, about six months of daily data).
#' @param n_starts Integer. Number of starting vectors (default 5). The
#'   shrinkage factors are spaced evenly from 1 to 0.45.
#' @param shrink Optional numeric vector in \eqn{(0, 1]} of shrinkage factors,
#'   overriding \code{n_starts} (expert setting; one start per factor).
#' @param stay Initial stay probability of each regime (default 0.95).
#'
#' @returns An object of class \code{"rsdc_starts"}: a list with
#' \describe{
#'   \item{starts}{\code{n_starts} \eqn{\times} \code{npar} matrix; each row is
#'     a feasible starting vector in the objective's packing order, with
#'     columns named as in \code{coef()}.}
#'   \item{loglik0}{Log-likelihood evaluated at each start (finite = feasible).}
#'   \item{method, N, K, p}{The model the starts are built for (validated by
#'     \code{\link{rsdc_estimate}} when the object is supplied).}
#'   \item{window, shrink, stay}{How the starts were constructed.}
#'   \item{regime_split}{Integer vector of length \eqn{T}: the quantile group
#'     (1 = lowest correlation) each observation was assigned to.}
#' }
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
#' st <- rsdc_starts(sim$observations, N = 2, method = "noX", n_starts = 3)
#' print(st)
#' fit <- rsdc_estimate("noX", residuals = sim$observations, N = 2,
#'                      control = list(start = st))
#' fit$start_logliks   # spread across starts (stability diagnostic)
#' }
#'
#' @seealso \code{\link{rsdc_estimate}} (accepts the object via
#'   \code{control$start}), \code{\link{rsdc_likelihood}} (parameter packing).
#' @importFrom stats cor quantile qlogis
#' @export
rsdc_starts <- function(residuals, N = 2, method = c("tvtp", "noX", "const"),
                        X = NULL, window = 126, n_starts = 5, shrink = NULL,
                        stay = 0.95) {
  method <- match.arg(method)
  if (!is.matrix(residuals) || ncol(residuals) < 2L)
    stop("residuals must be a numeric matrix with at least 2 columns.")
  if (any(!is.finite(residuals)))
    stop("residuals must not contain NA/NaN/Inf values.")
  Tn <- nrow(residuals); K <- ncol(residuals)
  if (method == "const") N <- 1L
  N <- as.integer(N)
  if (method != "const" && N < 2L) stop("N must be at least 2 for switching models.")
  if (method == "tvtp") {
    if (is.null(X)) stop("X must be provided for method = 'tvtp'.")
    X <- as.matrix(X)
    if (nrow(X) != Tn) stop("X must have the same number of rows as residuals.")
  }
  p <- if (method == "tvtp") ncol(X) else NA_integer_
  window <- as.integer(window)
  if (window < 20L) stop("window must be at least 20 observations.")
  if (Tn < max(2L * window, 30L * N))
    stop("Sample too short (T = ", Tn, ") for window = ", window,
         " and N = ", N, " regimes.")
  if (!is.numeric(stay) || stay <= 0 || stay >= 1)
    stop("stay must be in (0, 1).")
  if (is.null(shrink)) {
    n_starts <- as.integer(n_starts)
    if (n_starts < 1L) stop("n_starts must be at least 1.")
    shrink <- if (n_starts == 1L) 1 else seq(1, 0.45, length.out = n_starts)
  } else {
    if (any(shrink <= 0 | shrink > 1)) stop("shrink factors must be in (0, 1].")
    n_starts <- length(shrink)
  }

  ## 1. Rolling mean pairwise correlation, and quantile split into N groups
  lt <- lower.tri(diag(K))
  if (N == 1L) {
    regime_split <- rep(1L, Tn)
  } else {
    avg_corr <- rep(NA_real_, Tn)
    for (t in seq(window, Tn))
      avg_corr[t] <- mean(stats::cor(residuals[(t - window + 1L):t, ])[lt])
    avg_corr[seq_len(window - 1L)] <- avg_corr[window]
    qs <- stats::quantile(avg_corr, probs = seq(0, 1, length.out = N + 1L))
    if (anyDuplicated(qs))
      stop("Rolling correlations are too concentrated to split into ", N,
           " groups; reduce N or change window.")
    regime_split <- as.integer(cut(avg_corr, qs, labels = FALSE,
                                   include.lowest = TRUE))
  }

  ## 2. Empirical correlations per group (PD by construction), ascending order
  rho_groups <- lapply(seq_len(N), function(g)
    stats::cor(residuals[regime_split == g, , drop = FALSE])[lt])

  ## 3. Transition / covariate-coefficient initial values
  head_part <- switch(method,
    const = numeric(0),
    noX   = if (N == 2L) c(stay, stay) else {
      off <- (1 - stay) / (N - 1)
      unlist(lapply(seq_len(N), function(i) {
        row <- rep(off, N); row[i] <- stay
        row[seq_len(N - 1L)]                    # last column is the complement
      }))
    },
    tvtp  = {
      # Least-norm b with mean(X) %*% b = target logit, so the implied
      # transition matrix at the average covariate is persistent; covariate
      # slopes beyond that are left for the data to determine.
      xbar <- colMeans(X)
      if (sum(xbar^2) < 1e-10)
        stop("colMeans(X) is ~zero; include an intercept/constant column in X ",
             "(e.g. cbind(1, X)) so a persistent transition start can be built.")
      lam  <- xbar / sum(xbar^2)
      if (N == 2L) {
        rep(lam * stats::qlogis(stay), N)       # row-wise, identical rows
      } else {
        unlist(lapply(seq_len(N), function(i) {
          vecs <- matrix(0, N - 1L, p)          # softmax logits vs reference N
          if (i < N) vecs[i, ] <- lam * 3       # own state dominates
          else       vecs[]    <- rep(lam * -3, each = N - 1L)  # reference dominates
          as.vector(t(vecs))
        }))
      }
    })

  ## 4. One start per shrinkage factor (shrinking preserves PD)
  labs <- .rsdc_labels(method, N, K, p)
  starts <- t(vapply(shrink, function(s)
    c(head_part, unlist(rho_groups) * s), numeric(length(head_part) + N * K * (K - 1) / 2)))
  if (ncol(starts) != length(labs))
    stop("Internal error: start length (", ncol(starts),
         ") does not match the parameter labels (", length(labs), ").")
  colnames(starts) <- labs

  ## 5. Feasibility: log-likelihood at each start must be finite
  exog <- if (method == "tvtp") X else NULL
  loglik0 <- if (method == "const") {
    # Single-regime Gaussian log-likelihood (rsdc_likelihood needs N >= 2)
    apply(starts, 1, function(rho) {
      R <- diag(K); R[lower.tri(R)] <- rho; R[upper.tri(R)] <- t(R)[upper.tri(R)]
      ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
      if (any(ev <= 1e-10)) return(-Inf)
      sum(mvtnorm::dmvnorm(residuals, sigma = R, log = TRUE))
    })
  } else {
    apply(starts, 1, function(par)
      -rsdc_likelihood(par, y = residuals, exog = exog, K = K, N = N))
  }
  loglik0[loglik0 <= -1e9] <- -Inf

  structure(list(starts = starts, loglik0 = loglik0,
                 method = method, N = N, K = K, p = p,
                 window = window, shrink = shrink, stay = stay,
                 regime_split = regime_split, call = match.call()),
            class = "rsdc_starts")
}

#' @describeIn rsdc_starts Compact printer: model, dimensions, and feasibility
#'   of each start.
#' @param x An \code{"rsdc_starts"} object.
#' @param ... Unused; for generic compatibility.
#' @exportS3Method print rsdc_starts
print.rsdc_starts <- function(x, ...) {
  cat(sprintf("RSDC warm starts: method = \"%s\", N = %d, K = %d (npar = %d)\n",
              x$method, x$N, x$K, ncol(x$starts)))
  cat(sprintf("%d start(s), shrink = %s, window = %d\n",
              nrow(x$starts), paste(round(x$shrink, 2), collapse = ", "),
              x$window))
  feas <- is.finite(x$loglik0)
  cat(sprintf("Feasible: %d / %d | log-likelihood at starts: %s\n",
              sum(feas), length(feas),
              paste(round(x$loglik0, 1), collapse = " | ")))
  invisible(x)
}
