#' Inhomogeneous phase-type distributions
#'
#' Class of objects for inhomogeneous phase-type distributions.
#'
#' @slot name Name of the phase-type distribution.
#' @slot gfun A list comprising of the parameters.
#' @slot scale Scale.
#'
#' @return Class object.
#' @export
#'
setClass("iph",
  contains = c("ph"),
  slots = list(
    gfun = "list",
    scale = "numeric"
  )
)

#' Constructor function for inhomogeneous phase-type distributions
#'
#' @param ph An object of class \linkS4class{ph}.
#' @param alpha A probability vector.
#' @param S A sub-intensity matrix.
#' @param structure A valid ph structure.
#' @param dimension The dimension of the ph structure (if provided).
#' @param gfun Inhomogeneity transform.
#' @param gfun_pars The parameters of the inhomogeneity function.
#' @param scale Scale.
#'
#' @return An object of class \linkS4class{iph}.
#' @export
#'
#' @examples
#' iph(ph(structure = "coxian", dimension = 4), gfun = "pareto", gfun_pars = 3)
iph <- function(ph = NULL, gfun = NULL, gfun_pars = NULL, alpha = NULL, S = NULL, structure = NULL, dimension = 3, scale = 1) {
  if (all(is.null(c(gfun, gfun_pars)))) {
    stop("input inhomogeneity function and parameters")
  }
  if (is.null(ph)) {
    ph <- ph(alpha = alpha, S = S, structure = structure, dimension = dimension)
  }
  if (!gfun %in% c("pareto", "weibull", "lognormal", "loglogistic", "gompertz", "gev", "identity")) {
    stop("invalid gfun")
  }
  if (gfun %in% c("pareto", "weibull", "lognormal", "gompertz")) {
    if (is.null(gfun_pars)) gfun_pars <- 1
    if (length(gfun_pars) != 1 | sum(gfun_pars <= 0) > 0) {
      stop("gfun parameter should be positive and of length one")
    } else {
      names(gfun_pars) <- "beta"
    }
  }
  if (gfun %in% c("gev")) {
    if (is.null(gfun_pars)) gfun_pars <- c(0, 1, 1)
    if (length(gfun_pars) != 3 | (gfun_pars[2] > 0) == FALSE) {
      stop("gfun parameter should be of length three: mu, sigma, xi, and sigma > 0")
    } else {
      names(gfun_pars) <- c("mu", "sigma", "xi")
    }
  }
  if (gfun %in% c("loglogistic")) {
    if (is.null(gfun_pars)) gfun_pars <- c(1, 1)
    if (length(gfun_pars) != 2 | (gfun_pars[1] <= 0) | (gfun_pars[2] <= 0)) {
      stop("gfun parameter should be positive and of length two: alpha, theta > 0")
    } else {
      names(gfun_pars) <- c("alpha", "theta")
    }
  }
  f1 <- function(beta, t) t^(beta)
  f2 <- function(beta, t) log(t / beta + 1)
  f3 <- function(beta, t) log(t + 1)^(beta)
  f4 <- function(beta, t) log((t / beta[1])^(beta[2]) + 1)
  f5 <- function(beta, t) (exp(t * beta) - 1) / beta
  f6 <- function(beta, t, w) revers_data_trans(t, w, beta)
  nb <- which(gfun == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
  ginv <- base::eval(parse(text = paste("f", nb, sep = "")))

  f1 <- function(beta, t) t^(beta) * log(t)
  f2 <- function(beta, t) -t / (beta * t + beta^2)
  f3 <- function(beta, t) log(t + 1)^(beta) * log(log(t + 1))
  f4 <- NA
  f5 <- function(beta, t) exp(t * beta) * (t * beta - 1) / beta^2
  f6 <- NA
  nb <- which(gfun == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
  ginv_prime <- base::eval(parse(text = paste("f", nb, sep = "")))

  f1 <- function(beta, t) beta * t^(beta - 1)
  f2 <- function(beta, t) (t + beta)^(-1)
  f3 <- function(beta, t) beta * log(t + 1)^(beta - 1) / (t + 1)
  f4 <- NA
  f5 <- function(beta, t) exp(t * beta)
  f6 <- NA
  nb <- which(gfun == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
  lambda <- base::eval(parse(text = paste("f", nb, sep = "")))

  f1 <- function(beta, t) t^(beta - 1) + beta * t^(beta - 1) * log(t)
  f2 <- function(beta, t) -(t + beta)^(-2)
  f3 <- function(beta, t) log(t + 1)^(beta - 1) / (t + 1) + beta * log(t + 1)^(beta - 1) * log(log(t + 1)) / (t + 1)
  f4 <- NA
  f5 <- function(beta, t) t * exp(t * beta)
  f6 <- NA
  nb <- which(gfun == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
  lambda_prime <- base::eval(parse(text = paste("f", nb, sep = "")))

  name <- if (is(ph, "iph")) ph@name else paste("inhomogeneous ", ph@name, sep = "")

  methods::new("iph",
    name = name,
    pars = ph@pars,
    gfun = list(
      name = gfun, pars = gfun_pars, inverse = ginv,
      inverse_prime = ginv_prime, intensity = lambda, intensity_prime = lambda_prime
    ),
    scale = scale,
    fit = ph@fit
  )
}

#' Minimum method for inhomogeneous phase-type distributions
#'
#' @param x1 An object of class \linkS4class{iph}.
#' @param x2 An object of class \linkS4class{iph}.
#'
#' @return An object of class \linkS4class{iph}.
#' @export
#'
#' @examples
#' iph1 <- iph(ph(structure = "general", dimension = 3), gfun = "weibull", gfun_pars = 2)
#' iph2 <- iph(ph(structure = "gcoxian", dimension = 5), gfun = "weibull", gfun_pars = 2)
#' iph_min <- minimum(iph1, iph2)
#' iph_min
setMethod(
  "minimum", signature(x1 = "iph", x2 = "iph"),
  function(x1, x2) {
    if (x1@gfun$name != x2@gfun$name | !all(x1@gfun$pars == x2@gfun$pars)) {
      stop("g functions (or parameters) are different")
    }
    ph1 <- ph(alpha = x1@pars$alpha, S = x1@pars$S)
    ph2 <- ph(alpha = x2@pars$alpha, S = x2@pars$S)
    iph(ph = minimum(ph1, ph2), gfun = x1@gfun$name, gfun_pars = x1@gfun$pars)
  }
)

#' Maximum method for inhomogeneous phase-type distributions
#'
#' @param x1 An object of class \linkS4class{iph}.
#' @param x2 An object of class \linkS4class{iph}.
#'
#' @return An object of class \linkS4class{iph}.
#' @export
#'
#' @examples
#' iph1 <- iph(ph(structure = "general", dimension = 3), gfun = "weibull", gfun_pars = 2)
#' iph2 <- iph(ph(structure = "gcoxian", dimension = 5), gfun = "weibull", gfun_pars = 2)
#' iph_min <- maximum(iph1, iph2)
#' iph_min
setMethod(
  "maximum", signature(x1 = "iph", x2 = "iph"),
  function(x1, x2) {
    if (x1@gfun$name != x2@gfun$name | !all(x1@gfun$pars == x2@gfun$pars)) {
      stop("g functions (or parameters) are different")
    }
    ph1 <- ph(alpha = x1@pars$alpha, S = x1@pars$S)
    ph2 <- ph(alpha = x2@pars$alpha, S = x2@pars$S)
    iph(ph = maximum(ph1, ph2), gfun = x1@gfun$name, gfun_pars = x1@gfun$pars)
  }
)

#' Show method for inhomogeneous phase-type distributions
#'
#' @param object An object of class \linkS4class{iph}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "iph", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@pars)
  cat("g-function name: ", object@gfun$name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@gfun$pars)
})

#' Simulation method for inhomogeneous phase-type distributions
#'
#' @param x An object of class \linkS4class{iph}.
#' @param n An integer of length of realization.
#'
#' @return A realization of independent and identically distributed inhomogeneous
#'  phase-type variables.
#' @export
#'
#' @examples
#' obj <- iph(ph(structure = "general"), gfun = "lognormal", gfun_pars = 2)
#' sim(obj, n = 100)
setMethod("sim", c(x = "iph"), function(x, n = 1000) {
  name <- x@gfun$name
  pars <- x@gfun$pars
  scale <- x@scale
  if (name %in% c("pareto", "weibull", "lognormal", "loglogistic", "gompertz")) {
    U <- scale * riph(n, name, x@pars$alpha, x@pars$S, pars)
  }
  if (name %in% c("gev")) {
    U <- scale * rmatrixgev(n, x@pars$alpha, x@pars$S, pars[1], pars[2], pars[3])
  }
  U
})

#' Density method for inhomogeneous phase-type distributions
#'
#' @param x An object of class \linkS4class{iph}.
#' @param y A vector of locations.
#'
#' @return A vector containing the density evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- iph(ph(structure = "general"), gfun = "weibull", gfun_pars = 2)
#' dens(obj, c(1, 2, 3))
setMethod("dens", c(x = "iph"), function(x, y) {
  fn <- base::eval(parse(text = paste("m", x@gfun$name, "den", sep = "")))
  scale <- x@scale
  y_inf <- (y == Inf)
  dens <- y
  dens[!y_inf] <- fn(y / scale, x@pars$alpha, x@pars$S, x@gfun$pars) / scale
  dens[y_inf] <- 0
  dens
})

#' Distribution method for inhomogeneous phase-type distributions
#'
#' @param x An object of class \linkS4class{iph}.
#' @param q A vector of locations.
#' @param lower.tail Logical parameter specifying whether lower tail (CDF) or
#' upper tail is computed.
#'
#' @return A vector containing the CDF evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- iph(ph(structure = "general"), gfun = "weibull", gfun_pars = 2)
#' cdf(obj, c(1, 2, 3))
setMethod("cdf", c(x = "iph"), function(x,
                                        q,
                                        lower.tail = TRUE) {
  fn <- base::eval(parse(text = paste("m", x@gfun$name, "cdf", sep = "")))
  scale <- x@scale
  q_inf <- (q == Inf)
  cdf <- q
  cdf[!q_inf] <- fn(q[!q_inf] / scale, x@pars$alpha, x@pars$S, x@gfun$pars, lower.tail)
  cdf[q_inf] <- as.numeric(1 * lower.tail)
  cdf
})

#' Coef method for iph class
#'
#' @param object An object of class \linkS4class{iph}.
#'
#' @return Parameters of iph model.
#' @export
#'
#' @examples
#' obj <- iph(ph(structure = "general", dimension = 2), gfun = "lognormal", gfun_pars = 2)
#' coef(obj)
setMethod("coef", c(object = "iph"), function(object) {
  L <- object@pars
  L$gpars <- object@gfun$pars
  L
})
