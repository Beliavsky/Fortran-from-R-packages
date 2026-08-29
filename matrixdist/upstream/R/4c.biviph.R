#' Bivariate inhomogeneous phase-type distributions
#'
#' Class of objects for bivariate inhomogeneous phase-type distributions.
#'
#' @slot name Name of the phase type distribution.
#' @slot gfun A list comprising of the parameters.
#'
#' @return Class object.
#' @export
#'
setClass("biviph",
  contains = c("bivph"),
  slots = list(
    gfun = "list"
  )
)


#' Constructor function for bivariate inhomogeneous phase-type distributions
#'
#' @param bivph An object of class \linkS4class{bivph}.
#' @param alpha A probability vector.
#' @param S11 A sub-intensity matrix.
#' @param S12 A matrix.
#' @param S22 A sub-intensity matrix.
#' @param dimensions The dimensions of the bivariate phase-type (if no parameters are provided).
#' @param gfun Vector of inhomogeneity transforms.
#' @param gfun_pars List of parameters for the inhomogeneity functions.
#'
#' @return An object of class \linkS4class{biviph}.
#' @export
#'
#' @examples
#' under_bivph <- bivph(dimensions = c(3, 3))
#' biviph(under_bivph, gfun = c("weibull", "pareto"), gfun_pars = list(c(2), c(3)))
biviph <- function(bivph = NULL,
                   gfun = NULL,
                   gfun_pars = NULL,
                   alpha = NULL,
                   S11 = NULL,
                   S12 = NULL,
                   S22 = NULL,
                   dimensions = c(3, 3)) {
  if (all(is.null(c(gfun, gfun_pars)))) {
    stop("input inhomogeneity function and parameters")
  }
  d <- length(gfun)

  if (is.null(bivph)) {
    bivph <- bivph(alpha = alpha, S11 = S11, S12 = S12, S22 = S22, structure = structure, dimensions = dimensions)
  }

  if (!all(gfun %in% c("pareto", "weibull", "lognormal", "loglogistic", "gompertz", "gev", "identity"))) {
    stop("invalid gfun for at least one marginal")
  }

  if (all(gfun %in% c("pareto", "weibull", "lognormal", "gompertz"))) {
    for (i in 1:2) {
      if (is.null(gfun_pars[[i]])) gfun_pars[[i]] <- 1
      if (length(gfun_pars[[i]]) != 1 | sum(gfun_pars[[i]] <= 0) > 0) {
        stop(paste("gfun parameter for marginal", i, "should be positive and of length one"))
      } else {
        names(gfun_pars[[i]]) <- paste("beta", i, sep = "")
      }
    }
  }

  if (any(gfun %in% c("gev"))) {
    index_gev <- which(gfun == "gev")
    for (i in index_gev) {
      index_gev2 <- which(gfun[[i]] == "gev")
      if (is.null(gfun_pars[[i]])) gfun_pars[[i]] <- c(0, 1, 1)
      if (length(gfun_pars[[i]]) != 3 | (gfun_pars[[i]][2] > 0) == FALSE) {
        stop("gfun parameter should be of length three: mu, sigma, xi, and sigma > 0")
      } else {
        names(gfun_pars) <- c("mu", "sigma", "xi")
      }
    }
  }

  if (any(gfun == "loglogistic")) {
    index_log <- which(gfun == "loglogistic")
    for (i in index_log) {
      if (is.null(gfun_pars[[i]])) gfun_pars[[i]] <- c(1, 1)
      if (length(gfun_pars[[i]]) != 2 | (gfun_pars[[i]][1] <= 0) | (gfun_pars[[i]][2] <= 0)) {
        stop("gfun parameter should be positive and of length two: alpha, theta > 0")
      } else {
        names(gfun_pars[[i]]) <- c(paste("alpha", i, sep = ""), paste("theta", i, sep = ""))
      }
    }
  }

  ginv <- list()
  ginv_prime <- list()
  lambda <- list()
  lambda_prime <- list()

  for (i in 1:2) {
    f1 <- function(beta, t) t^(beta)
    f2 <- function(beta, t) log(t / beta + 1)
    f3 <- function(beta, t) log(t + 1)^(beta)
    f4 <- function(beta, t) log((t / beta[1])^(beta[2]) + 1)
    f5 <- function(beta, t) (exp(t * beta) - 1) / beta
    f6 <- function(beta, t, w) revers_data_trans(t, w, beta)
    nb <- which(gfun[i] == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
    ginv[[i]] <- base::eval(parse(text = paste("f", nb, sep = "")))

    f1 <- function(beta, t) t^(beta) * log(t)
    f2 <- function(beta, t) -t / (beta * t + beta^2)
    f3 <- function(beta, t) log(t + 1)^(beta) * log(log(t + 1))
    f4 <- NA
    f5 <- function(beta, t) exp(t * beta) * (t * beta - 1) / beta^2
    f6 <- NA
    nb <- which(gfun[i] == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
    ginv_prime[[i]] <- base::eval(parse(text = paste("f", nb, sep = "")))

    f1 <- function(beta, t) beta * t^(beta - 1)
    f2 <- function(beta, t) (t + beta)^(-1)
    f3 <- function(beta, t) beta * log(t + 1)^(beta - 1) / (t + 1)
    f4 <- NA
    f5 <- function(beta, t) exp(t * beta)
    f6 <- NA
    nb <- which(gfun[i] == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
    lambda[[i]] <- base::eval(parse(text = paste("f", nb, sep = "")))

    f1 <- function(beta, t) t^(beta - 1) + beta * t^(beta - 1) * log(t)
    f2 <- function(beta, t) -(t + beta)^(-2)
    f3 <- function(beta, t) log(t + 1)^(beta - 1) / (t + 1) + beta * log(t + 1)^(beta - 1) * log(log(t + 1)) / (t + 1)
    f4 <- NA
    f5 <- function(beta, t) t * exp(t * beta)
    f6 <- NA
    nb <- which(gfun[i] == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
    lambda_prime[[i]] <- base::eval(parse(text = paste("f", nb, sep = "")))
  }
  name <- if (is(bivph, "biviph")) bivph@name else paste("inhomogeneous ", bivph@name, sep = "")

  methods::new("biviph",
    name = name,
    pars = bivph@pars,
    gfun = list(
      name = gfun,
      pars = gfun_pars,
      inverse = ginv,
      inverse_prime = ginv_prime,
      intensity = lambda,
      intensity_prime = lambda_prime
    )
  )
}

#' Show method for bivariate inhomogeneous phase-type distributions
#'
#' @param object An object of class \linkS4class{biviph}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "biviph", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@pars)
  cat("g-function name:", object@gfun$name, "\n")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@gfun$pars)
})

#' Simulation method for bivariate inhomogeneous phase-type distributions
#'
#' @param x An object of class \linkS4class{biviph}.
#' @param n An integer of length of realization.
#'
#' @return A realization of independent and identically distributed bivariate
#'  inhomogeneous phase-type vector.
#' @export
#'
#' @examples
#' under_bivph <- bivph(dimensions = c(3, 3))
#' obj <- biviph(under_bivph, gfun = c("weibull", "pareto"), gfun_pars = list(c(2), c(3)))
#' sim(obj, n = 100)
setMethod("sim", c(x = "biviph"), function(x, n = 1000) {
  p1_aux <- dim(x@pars$S11)[1]
  p2_aux <- dim(x@pars$S22)[1]
  alpha_aux <- c(x@pars$alpha, rep(0, p2_aux))
  S_aux <- merge_matrices(x@pars$S11, x@pars$S12, x@pars$S22)
  R_aux <- matrix(c(c(rep(1, p1_aux), rep(0, p2_aux)), c(rep(0, p1_aux), rep(1, p2_aux))), ncol = 2)
  rMIPHstar(n, alpha_aux, S_aux, R_aux, x@gfun$name, x@gfun$pars)
})

#' Marginal method for biviph class
#'
#' @param x An object of class \linkS4class{biviph}.
#' @param mar Indicator of which marginal.
#' @return An object of the of class \linkS4class{iph}.
#' @export
#'
#' @examples
#' under_bivph <- bivph(dimensions = c(3, 3))
#' obj <- biviph(under_bivph, gfun = c("weibull", "pareto"), gfun_pars = list(c(2), c(3)))
#' marginal(obj, 1)
setMethod("marginal", c(x = "biviph"), function(x, mar = 1) {
  if (!(mar %in% 1:2)) {
    stop("maringal provided not available")
  }
  if (mar == 1) {
    x0 <- iph(ph(alpha = x@pars$alpha, S = x@pars$S11), gfun = x@gfun$name[mar], gfun_pars = x@gfun$pars[[mar]])
  } else {
    alpha0 <- x@pars$alpha %*% base::solve(-x@pars$S11) %*% x@pars$S12
    x0 <- iph(ph(alpha = alpha0, S = x@pars$S22), gfun = x@gfun$name[mar], gfun_pars = x@gfun$pars[[mar]])
  }
  x0
})

#' Density method for bivariate inhomogeneous phase-type distributions
#'
#' @param x An object of class \linkS4class{biviph}.
#' @param y A matrix of locations.
#'
#' @return A vector containing the joint density evaluations at the given locations.
#' @export
#'
#' @examples
#' under_bivph <- bivph(dimensions = c(3, 3))
#' obj <- biviph(under_bivph, gfun = c("weibull", "pareto"), gfun_pars = list(c(2), c(3)))
#' dens(obj, matrix(c(0.5, 1), ncol = 2))
setMethod("dens", c(x = "biviph"), function(x, y) {
  if (is.vector(y)) {
    y <- t(y)
  }
  gfun.pars <- x@gfun$pars
  y1_inv <- x@gfun$inverse[[1]](gfun.pars[[1]], y[, 1])
  y1_int <- x@gfun$intensity[[1]](gfun.pars[[1]], y[, 1])
  y2_inv <- x@gfun$inverse[[2]](gfun.pars[[2]], y[, 2])
  y2_int <- x@gfun$intensity[[2]](gfun.pars[[2]], y[, 2])
  dens <- bivph_density(cbind(y1_inv, y2_inv), x@pars$alpha, x@pars$S11, x@pars$S12, x@pars$S22) * y1_int * y2_int
  unname(dens)
})

#' Coef method for biviph class
#'
#' @param object An object of class \linkS4class{biviph}.
#'
#' @return Parameters of bivariate inhomogeneous phase-type model.
#' @export
#'
#' @examples
#' under_bivph <- bivph(dimensions = c(3, 3))
#' obj <- biviph(under_bivph, gfun = c("weibull", "pareto"), gfun_pars = list(c(2), c(3)))
#' coef(obj)
setMethod("coef", c(object = "biviph"), function(object) {
  L <- object@pars
  L$gfun <- object@gfun$name
  L$gpars <- object@gfun$pars
  L
})