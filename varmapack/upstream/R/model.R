#' Create a VARMA or VARMAX model
#'
#' Create a model object for simulation of a Gaussian VAR, VMA, VARMA, or
#' VARMAX time series. The object stores model parameters and provides a
#' `$sim()` method for generating one or more independent series.
#'
#' @param A Autoregressive coefficient matrices as an `r` by `r` matrix or an
#'   `r` by `r` by `p` array. `NULL` gives no autoregressive terms.
#' @param B Moving-average coefficient matrices as an `r` by `r` matrix or an
#'   `r` by `r` by `q` array. `NULL` gives no moving-average terms.
#' @param C Exogenous coefficient matrices as an `r` by `d` matrix or an
#'   `r` by `d` by `s` array. `NULL` gives a VARMA model rather than VARMAX.
#' @param Sig `r` by `r` innovation covariance matrix.
#' @param mu Optional `r` vector or `r` by `nmu` matrix of time-series means.
#'   The last mean vector repeats for the rest of each series. Means are not
#'   supported for VARMAX models.
#'
#' @return A `VarmapackModel` object.
#'
#' @section Simulation:
#' `model$sim(length, nrep = 1L, X0 = NULL, z = NULL, rng = NULL,
#' return_shocks = FALSE)` simulates `nrep` independent series of the supplied
#' `length`.
#'
#' For VARMA models, `X0` is optional. If provided, it is an `r` by `nX0`
#' matrix or an `r` by `nX0` by `nrep` array. Its second dimension must be at
#' least `max(p, q)`, and its third dimension, when present, must be 1 or
#' `nrep`.
#' Nonstationary VARMA models require `X0`; with MA terms, `Sig` must be
#' positive definite and startup shocks are conditioned on the residual
#' equations implied by the supplied history.
#'
#' For VARMAX models, `z` is required. `X0` has the same form and its second
#' dimension must be at least `max(p, q, s - 1)`. It may be omitted when this
#' minimum is zero. The exogenous input `z` is a `d` by `length` matrix or a `d`
#' by `length` by `nrep` array. Its third dimension, when present, must be 1 or
#' `nrep`.
#'
#' Pass a `randompack::randompack_rng()` object through `rng` to control the
#' random stream. When it is omitted, Varmapack uses a temporary randomized
#' default Randompack generator. With `return_shocks = FALSE`, the result is an
#' `r` by `length` by `nrep` array. With `return_shocks = TRUE`, it is a list
#' with components `X` and `E`, each with that shape.
#'
#' @section Model analysis:
#' `model$acvf(maxlag)` returns theoretical VARMA autocovariances,
#' `model$psi(maxlag)` returns impulse-response matrices, and
#' `model$irf(maxlag)` returns orthogonalized impulse-response matrices.
#' `model$specrad()` and `model$ma_specrad()` return the AR and MA spectral
#' radii, respectively. Theoretical autocovariances are not available for
#' VARMAX models with exogenous terms.
#'
#' @examples
#' A <- matrix(c(0.4, 0.1, -0.2, 0.3), 2, 2)
#' model <- varmapack_model(A = A, Sig = diag(2))
#' X <- model$sim(100)
#'
#' rng <- randompack::randompack_rng()
#' rng$seed(123)
#' X <- model$sim(100, nrep = 3, rng = rng)
#'
#' @export
varmapack_model <- function(A = NULL, B = NULL, C = NULL, Sig, mu = NULL) {
  VarmapackModel$new(A = A, B = B, C = C, Sig = Sig, mu = mu)
}

.vp_numeric_array <- function(x, name) {
  if (!is.numeric(x) || anyNA(x) || !all(is.finite(x)))
    stop(name, " must be a finite numeric array", call. = FALSE)
  x
}

.vp_square_terms <- function(x, name, r) {
  if (is.null(x))
    return(NULL)
  x <- .vp_numeric_array(x, name)
  d <- dim(x)
  if (is.null(d) || !(length(d) %in% c(2L, 3L)))
    stop(name, " must be a matrix or three-dimensional array", call. = FALSE)
  if (d[1] != r || d[2] != r || any(d <= 0L))
    stop(name, " must have dimensions r by r by p", call. = FALSE)
  if (length(d) == 2L)
    d <- c(d, 1L)
  array(as.double(x), dim = d)
}

.vp_exogenous_terms <- function(C, r) {
  if (is.null(C))
    return(NULL)
  C <- .vp_numeric_array(C, "C")
  d <- dim(C)
  if (is.null(d) || !(length(d) %in% c(2L, 3L)))
    stop("C must be a matrix or three-dimensional array", call. = FALSE)
  if (d[1] != r || any(d <= 0L))
    stop("C must have dimensions r by d by s", call. = FALSE)
  if (length(d) == 2L)
    d <- c(d, 1L)
  array(as.double(C), dim = d)
}

.vp_mean_path <- function(mu, r) {
  if (is.null(mu))
    return(NULL)
  mu <- .vp_numeric_array(mu, "mu")
  d <- dim(mu)
  if (is.null(d)) {
    if (length(mu) != r)
      stop("mu must have length r or dimensions r by nmu", call. = FALSE)
    return(matrix(as.double(mu), nrow = r))
  }
  if (length(d) != 2L || d[1] != r || d[2] <= 0L)
    stop("mu must have length r or dimensions r by nmu", call. = FALSE)
  matrix(as.double(mu), nrow = r, ncol = d[2])
}

.vp_count <- function(x, name, minimum) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max)
    stop(name, " must be an integer at least ", minimum, call. = FALSE)
  as.integer(x)
}

.vp_start_values <- function(X0, r, nrep, minimum, required) {
  if (is.null(X0)) {
    if (required)
      stop("X0 is required for VARMAX simulation", call. = FALSE)
    return(NULL)
  }
  X0 <- .vp_numeric_array(X0, "X0")
  d <- dim(X0)
  if (is.null(d) || !(length(d) %in% c(2L, 3L)))
    stop("X0 must be a matrix or three-dimensional array", call. = FALSE)
  if (d[1] != r || d[2] < minimum || any(d <= 0L))
    stop("X0 has invalid dimensions", call. = FALSE)
  npaths <- if (length(d) == 2L) 1L else d[3]
  if (!(npaths %in% c(1L, nrep)))
    stop("X0 must contain one or nrep startup paths", call. = FALSE)
  array(as.double(X0), dim = d)
}

.vp_exogenous_input <- function(z, d, n, nrep) {
  if (is.null(z))
    stop("z is required for VARMAX simulation", call. = FALSE)
  z <- .vp_numeric_array(z, "z")
  dz <- dim(z)
  if (is.null(dz) || !(length(dz) %in% c(2L, 3L)))
    stop("z must be a matrix or three-dimensional array", call. = FALSE)
  if (dz[1] != d || dz[2] != n || any(dz <= 0L))
    stop("z has invalid dimensions", call. = FALSE)
  npaths <- if (length(dz) == 2L) 1L else dz[3]
  if (!(npaths %in% c(1L, nrep)))
    stop("z must contain one or nrep input paths", call. = FALSE)
  array(as.double(z), dim = dz)
}

VarmapackModel <- R6::R6Class(
  "VarmapackModel",
  public = list(
    A = NULL,
    B = NULL,
    C = NULL,
    Sig = NULL,
    mu = NULL,
    p = NULL,
    q = NULL,
    r = NULL,
    s = NULL,
    d = NULL,
    initialize = function(A = NULL, B = NULL, C = NULL, Sig, mu = NULL) {
      Sig <- .vp_numeric_array(Sig, "Sig")
      if (!is.matrix(Sig) || nrow(Sig) <= 0L || nrow(Sig) != ncol(Sig))
        stop("Sig must be a non-empty square matrix", call. = FALSE)
      if (!isTRUE(all(Sig == t(Sig))))
        stop("Sig must be symmetric", call. = FALSE)
      self$r <- nrow(Sig)
      self$A <- .vp_square_terms(A, "A", self$r)
      self$B <- .vp_square_terms(B, "B", self$r)
      self$C <- .vp_exogenous_terms(C, self$r)
      self$Sig <- matrix(as.double(Sig), nrow = self$r, ncol = self$r)
      self$mu <- .vp_mean_path(mu, self$r)
      if (!is.null(self$C) && !is.null(self$mu))
        stop("mu is not supported for VARMAX models", call. = FALSE)
      self$p <- if (is.null(self$A)) 0L else dim(self$A)[3]
      self$q <- if (is.null(self$B)) 0L else dim(self$B)[3]
      self$s <- if (is.null(self$C)) 0L else dim(self$C)[3]
      self$d <- if (is.null(self$C)) 0L else dim(self$C)[2]
    },
    sim = function(length, nrep = 1L, X0 = NULL, z = NULL, rng = NULL,
                   return_shocks = FALSE) {
      length <- .vp_count(length, "length", 1L)
      nrep <- .vp_count(nrep, "nrep", 1L)
      if (!is.logical(return_shocks) || length(return_shocks) != 1L ||
          is.na(return_shocks))
        stop("return_shocks must be TRUE or FALSE", call. = FALSE)
      if (!is.null(rng) && (!inherits(rng, "RandompackRNG") || is.null(rng$ptr)))
        stop("rng must be a randompack_rng() object", call. = FALSE)
      if (is.null(self$C)) {
        if (!is.null(z))
          stop("z can only be supplied for VARMAX models", call. = FALSE)
        X0 <- .vp_start_values(X0, self$r, nrep, max(self$p, self$q), FALSE)
      }
      else {
        minimum <- max(self$p, self$q, self$s - 1L)
        X0 <- .vp_start_values(X0, self$r, nrep, minimum, minimum > 0L)
        z <- .vp_exogenous_input(z, self$d, length, nrep)
      }
      .Call("varmapack_sim_R", self$A, self$B, self$C, self$Sig, self$mu, X0, z,
            length, nrep, return_shocks, if (is.null(rng)) NULL else rng$ptr,
            PACKAGE = "varmapack")
    },
    acvf = function(maxlag) {
      if (!is.null(self$C))
        stop("theoretical autocovariances are available only for VARMA models",
             call. = FALSE)
      maxlag <- .vp_count(maxlag, "maxlag", 0L)
      .Call("varmapack_acvf_R", self$A, self$B, self$Sig, maxlag,
            PACKAGE = "varmapack")
    },
    psi = function(maxlag) {
      maxlag <- .vp_count(maxlag, "maxlag", 0L)
      .Call("varmapack_psi_R", self$A, self$B, self$Sig, maxlag,
            PACKAGE = "varmapack")
    },
    irf = function(maxlag) {
      maxlag <- .vp_count(maxlag, "maxlag", 0L)
      .Call("varmapack_irf_R", self$A, self$B, self$Sig, maxlag,
            PACKAGE = "varmapack")
    },
    specrad = function() {
      .Call("varmapack_specrad_R", self$A, self$Sig, PACKAGE = "varmapack")
    },
    ma_specrad = function() {
      .Call("varmapack_ma_specrad_R", self$B, self$Sig, PACKAGE = "varmapack")
    }
  )
)
