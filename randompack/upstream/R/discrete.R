# R/randompack_discrete.R

methods_discrete <- list(
  int = function(len, min, max) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    .Call("randompack_int_R", self$ptr, len, min, max,
      PACKAGE = "randompack")
  },

  perm = function(n) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    .Call("randompack_perm_R", self$ptr, n, PACKAGE = "randompack")
  },

  sample = function(n, k) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    .Call("randompack_sample_R", self$ptr, n, k, PACKAGE = "randompack")
  },

  raw = function(n) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    .Call("randompack_raw_R", self$ptr, n, PACKAGE = "randompack")
  }
)

NULL

NULL
