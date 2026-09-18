#' Sample autocovariances of an observed series
#'
#' Compute sample autocovariance matrices up to a specified lag for a numeric
#' time-series matrix with variables in rows and observations in columns.
#'
#' @param X Numeric `r` by `n` observed time-series matrix.
#' @param maxlag Largest lag to compute, between zero and `n - 1`.
#' @param norm Either `"ML"`, which divides every lag by `n`, or `"C"`, which
#'   divides the lag-`k` covariance by `n - k`.
#'
#' @return An `r` by `r` by `maxlag + 1` array. Its `k + 1` plane estimates
#'   `Cov(x_t, x_{t-k})`.
#'
#' @examples
#' X <- rbind(1:10, (1:10)^2)
#' varmapack_autocov(X, maxlag = 2)
#'
#' @export
varmapack_autocov <- function(X, maxlag, norm = "ML") {
  X <- .vp_numeric_array(X, "X")
  if (!is.matrix(X) || nrow(X) <= 0L || ncol(X) <= 0L)
    stop("X must be a non-empty numeric matrix", call. = FALSE)
  maxlag <- .vp_count(maxlag, "maxlag", 0L)
  if (maxlag >= ncol(X))
    stop("maxlag must be less than the number of observations", call. = FALSE)
  if (!is.character(norm) || length(norm) != 1L || is.na(norm) ||
      !(toupper(norm) %in% c("ML", "C")))
    stop("norm must be \"ML\" or \"C\"", call. = FALSE)
  .Call("varmapack_autocov_R", matrix(as.double(X), nrow = nrow(X)), maxlag,
        toupper(norm), PACKAGE = "varmapack")
}

#' Convert covariances to correlations
#'
#' Convert a covariance matrix or sequence of autocovariance matrices to
#' correlations using the marginal standard deviations at lag zero.
#'
#' @param cov A finite numeric `r` by `r` matrix or `r` by `r` by
#'   `maxlag + 1` array. The diagonal entries of its lag-zero matrix must be
#'   positive.
#'
#' @return A numeric array with the same dimensions as `cov`. Lag-zero
#'   diagonal entries are exactly one. Other entries are not clipped to
#'   `[-1,1]`.
#'
#' @examples
#' cov <- array(c(4, 3, 3, 9, 2, 6, -3, 1.5), dim = c(2, 2, 2))
#' varmapack_cov2corr(cov)
#'
#' @export
varmapack_cov2corr <- function(cov) {
  cov <- .vp_numeric_array(cov, "cov")
  d <- dim(cov)
  if (is.null(d) || !(length(d) %in% c(2L, 3L)))
    stop("cov must be a matrix or three-dimensional array", call. = FALSE)
  if (d[1] <= 0L || d[2] != d[1] || any(d <= 0L))
    stop("cov must have dimensions r by r by maxlag + 1", call. = FALSE)
  cov <- array(as.double(cov), dim = d)
  .Call("varmapack_cov2corr_R", cov, PACKAGE = "varmapack")
}

#' Create a VARMA testcase model
#'
#' Create a built-in named testcase or an unnamed random, deterministic, or
#' spectral-radius-controlled VARMA testcase.
#'
#' @param which A built-in testcase name, a one-based built-in testcase index,
#'   or one of `"random"`, `"deterministic"`, and `"rho"`.
#' @param p,q,r Required AR order, MA order, and series dimension for unnamed
#'   testcases.
#' @param rho Target spectral radius when `which = "rho"`.
#' @param rng Optional `randompack::randompack_rng()` used for a random testcase.
#'
#' @return A `VarmapackModel` object.
#'
#' @examples
#' model <- varmapack_testcase("smallARMA1")
#' model <- varmapack_testcase("rho", p = 3, q = 1, r = 2, rho = 0.8)
#'
#' @export
varmapack_testcase <- function(which = "random", p = NULL, q = NULL, r = NULL,
                               rho = 0, rng = NULL) {
  index <- 0L
  name <- ""
  unnamed <- FALSE
  if (is.character(which) && length(which) == 1L && !is.na(which)) {
    if (which == "random") {
      unnamed <- TRUE
    }
    else if (which == "deterministic") {
      index <- -1L
      unnamed <- TRUE
    }
    else if (which == "rho") {
      name <- "rho"
      unnamed <- TRUE
    }
    else if (which == "max") {
      stop("use varmapack_testcases() to list named testcases", call. = FALSE)
    }
    else {
      name <- which
    }
  }
  else if (is.numeric(which) && length(which) == 1L && !is.na(which) &&
           is.finite(which) && which == floor(which) && which >= 1L &&
           which <= .Machine$integer.max) {
    index <- as.integer(which)
  }
  else {
    stop("which must be a testcase name or positive index", call. = FALSE)
  }
  if (unnamed) {
    p <- .vp_count(p, "p", 0L)
    q <- .vp_count(q, "q", 0L)
    r <- .vp_count(r, "r", 1L)
  }
  else {
    p <- q <- r <- 0L
  }
  if (!is.numeric(rho) || length(rho) != 1L || is.na(rho) || !is.finite(rho) || rho < 0)
    stop("rho must be a finite non-negative number", call. = FALSE)
  if (!is.null(rng) && (!inherits(rng, "RandompackRNG") || is.null(rng$ptr)))
    stop("rng must be a randompack_rng() object", call. = FALSE)
  parts <- .Call("varmapack_testcase_R", name, index, as.double(rho), p, q, r,
                 if (is.null(rng)) NULL else rng$ptr, PACKAGE = "varmapack")
  varmapack_model(A = parts$A, B = parts$B, Sig = parts$Sig)
}

#' List named VARMA testcases
#'
#' @return A data frame with columns `index`, `name`, `p`, `q`, and `r`.
#'
#' @examples
#' varmapack_testcases()
#'
#' @export
varmapack_testcases <- function() {
  parts <- .Call("varmapack_testcases_R", PACKAGE = "varmapack")
  data.frame(parts, stringsAsFactors = FALSE)
}
