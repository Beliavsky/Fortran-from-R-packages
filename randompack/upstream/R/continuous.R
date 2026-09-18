# R/randompack_continuous.R

methods_continuous <- list(
  unif = function(len, a = 0, b = 1) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    .Call("randompack_unif_R", self$ptr, len, a, b, PACKAGE = "randompack")
  },

  normal = function(len, mu = 0, sigma = 1) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    .Call("randompack_normal_R", self$ptr, len, mu, sigma,
          PACKAGE = "randompack")
  },

  skew_normal = function(len, mu = 0, sigma = 1, alpha) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    .Call("randompack_skew_normal_R", self$ptr, len, mu, sigma, alpha,
          PACKAGE = "randompack")
  },

  lognormal = function(len, mu = 0, sigma = 1) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    .Call("randompack_lognormal_R", self$ptr, len, mu, sigma,
          PACKAGE = "randompack")
  },

  gumbel = function(len, mu = 0, beta = 1) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    .Call("randompack_gumbel_R", self$ptr, len, mu, beta,
          PACKAGE = "randompack")
  },

  pareto = function(len, xm, alpha) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    .Call("randompack_pareto_R", self$ptr, len, xm, alpha,
          PACKAGE = "randompack")
  },

  exp = function(len, scale = 1) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    .Call("randompack_exp_R", self$ptr, len, scale, PACKAGE = "randompack")
  },

  gamma = function(len, shape, scale = 1) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    .Call("randompack_gamma_R", self$ptr, len, shape, scale,
          PACKAGE = "randompack")
  },

  chi2 = function(len, nu) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    .Call("randompack_chi2_R", self$ptr, len, nu, PACKAGE = "randompack")
  },

  beta = function(len, a, b) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    .Call("randompack_beta_R", self$ptr, len, a, b, PACKAGE = "randompack")
  },

  t = function(len, nu) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    .Call("randompack_t_R", self$ptr, len, nu, PACKAGE = "randompack")
  },

  f = function(len, nu1, nu2) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    .Call("randompack_f_R", self$ptr, len, nu1, nu2, PACKAGE = "randompack")
  },

  weibull = function(len, shape, scale = 1) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    .Call("randompack_weibull_R", self$ptr, len, shape, scale,
          PACKAGE = "randompack")
  },
  mvn = function(n, Sigma, mu = NULL) {
    if (is.null(self$ptr)) stop("RNG is not initialized")
    if (!is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n))
      stop("n must be a finite numeric scalar")
    if (n < 0 || n != floor(n)) stop("n must be a non-negative integer")
    n <- as.integer(n)
    if (is.null(Sigma)) stop("Sigma must not be NULL")
    if (!is.matrix(Sigma) || !is.numeric(Sigma))
      stop("Sigma must be a numeric matrix")
    d <- nrow(Sigma)
    if (d <= 0L || ncol(Sigma) != d)
      stop("Sigma must be a non-empty square matrix")
    if (anyNA(Sigma)) stop("Sigma must not contain NA")
    if (!all(is.finite(Sigma))) stop("Sigma must be finite")
    if (!all(Sigma == t(Sigma))) stop("Sigma must be symmetric")
    storage.mode(Sigma) <- "double"
    if (!is.null(mu)) {
      if (!is.numeric(mu)) stop("mu must be a numeric vector")
      mu <- as.numeric(mu)
      if (length(mu) != d) stop("mu length must match Sigma")
      if (anyNA(mu)) stop("mu must not contain NA")
      if (!all(is.finite(mu))) stop("mu must be finite")
    }
    .Call("randompack_mvn_R", self$ptr, n, Sigma, mu,
          PACKAGE = "randompack")
  }
)
