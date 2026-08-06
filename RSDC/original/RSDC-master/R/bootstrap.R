#' Parametric bootstrap standard errors for a fitted RSDC model
#'
#' Simulates \code{B} data sets from the fitted data-generating process, re-estimates
#' each one — \emph{warm-started} at the original MLE so the global search is skipped —
#' and returns the bootstrap covariance, standard errors and percentile intervals of
#' the parameters. This complements the analytic (Hessian/OPG/sandwich) covariances
#' from \code{\link{vcov.rsdc_fit}} and, being simulation-based, respects the parameter
#' bounds (\eqn{\rho \in (-1,1)}, transition probabilities in \eqn{[0,1]}).
#'
#' @details
#' For \code{method = "tvtp"} the data are simulated with \code{\link{rsdc_simulate}}
#' from the fitted \code{beta} and covariate matrix \code{X}; for \code{"noX"} from the
#' fitted (fixed) transition matrix; for \code{"const"} as i.i.d. draws from the single
#' regime. Each replicate is re-estimated with the same \code{method}/\code{N} and
#' \code{control = list(start = coef(object), compute_se = FALSE)}, which warm-starts the
#' local optimizer at the original estimates (fast, deterministic given a seed).
#' Replicates whose re-estimation fails are dropped.
#'
#' @param object An object of class \code{"rsdc_fit"}.
#' @param B Integer number of bootstrap replications (default 199).
#' @param X Covariate matrix for \code{"tvtp"} (defaults to the stored estimation covariates).
#' @param seed Optional integer seed for reproducibility.
#' @param level Confidence level for the percentile intervals (default 0.95).
#' @param cores Integer (default 1). Number of cores used to re-estimate the
#'   replicates in parallel (via the \pkg{parallel} package). All \code{B}
#'   data sets are simulated first, from a single RNG stream, so the result
#'   is identical for any number of cores.
#'
#' @return A list with \code{replicates} (\eqn{B' \times} npar matrix of successful
#'   re-estimates), \code{vcov}, \code{se}, \code{ci} (percentile intervals) and the
#'   effective number of replicates \code{B}.
#'
#' @seealso \code{\link{vcov.rsdc_fit}}, \code{\link{rsdc_simulate}}, \code{\link{rsdc_estimate}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- cbind(1, as.numeric(scale(seq_len(300))))
#' y <- scale(matrix(rnorm(300 * 2), 300, 2))
#' fit <- rsdc_estimate("tvtp", residuals = y, N = 2, X = X)
#' bs <- rsdc_bootstrap(fit, B = 50, seed = 1)
#' bs$se
#' bs$ci
#' }
#'
#' @importFrom stats cov sd quantile complete.cases
#' @export
rsdc_bootstrap <- function(object, B = 199, X = NULL, seed = NULL, level = 0.95,
                           cores = 1) {
  stopifnot(inherits(object, "rsdc_fit"))
  if (is.null(object$par)) stop("object has no parameter vector to bootstrap.")
  if (!is.null(seed)) set.seed(seed)

  method <- object$method; N <- object$N
  if (method == "tvtp") {
    if (is.null(X)) X <- object$X
    if (is.null(X)) stop("Provide X for a tvtp fit (covariates were not stored).")
  }

  labs <- names(object$coefficients)
  ctrl <- list(start = object$par, compute_se = FALSE)

  # Simulate all B data sets first, from one uninterrupted RNG stream: the
  # warm-started re-estimations below are deterministic, so they can run in
  # any order (and in parallel) without affecting the result.
  # tvtp uses nrow(X); noX/const default to the fitted length (object$nobs).
  sims <- lapply(seq_len(B), function(b)
    tryCatch(simulate(object, X = X), error = function(e) NULL))

  refit_one <- function(sim) {
    if (is.null(sim)) return(NULL)
    ysim <- scale(sim$observations)                       # unit-variance residuals
    fitb <- tryCatch(
      switch(method,
             tvtp  = rsdc_estimate("tvtp", residuals = ysim, N = N, X = X, control = ctrl),
             noX   = rsdc_estimate("noX",  residuals = ysim, N = N,        control = ctrl),
             const = rsdc_estimate("const", residuals = ysim,             control = ctrl)),
      error = function(e) NULL)
    if (!is.null(fitb) && length(fitb$par) == length(labs)) fitb$par else NULL
  }
  pars <- .rsdc_lapply(sims, refit_one, cores = cores)

  ok   <- !vapply(pars, is.null, logical(1))
  reps <- do.call(rbind, pars[ok])
  if (is.null(reps)) reps <- matrix(NA_real_, 0, length(labs))
  colnames(reps) <- labs
  reps <- reps[stats::complete.cases(reps), , drop = FALSE]
  if (nrow(reps) < 2L)
    stop("Too few successful bootstrap replicates (", nrow(reps), "); increase B.")

  V  <- stats::cov(reps)
  se <- apply(reps, 2, stats::sd)
  a  <- (1 - level) / 2
  ci <- t(apply(reps, 2, stats::quantile, probs = c(a, 1 - a), na.rm = TRUE))
  dimnames(V) <- list(labs, labs); names(se) <- labs; rownames(ci) <- labs
  list(replicates = reps, vcov = V, se = se, ci = ci, B = nrow(reps))
}
