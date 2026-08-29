#' Multivariate inhomogeneous phase-type distributions
#'
#' Class of objects for multivariate inhomogeneous phase-type distributions.
#'
#' @slot name Name of the phase type distribution.
#' @slot gfun A list comprising of the parameters.
#' @slot scale Scale.
#'
#' @return Class object.
#' @export
#'
setClass("miph",
  contains = c("mph"),
  slots = list(
    gfun = "list",
    scale = "numeric"
  )
)

#' Constructor function for multivariate inhomogeneous phase-type distributions
#'
#' @param mph An object of class \linkS4class{mph}.
#' @param alpha A probability vector.
#' @param S A list of sub-intensity matrices.
#' @param structure A vector of valid ph structures.
#' @param dimension The dimension of the ph structure (if provided).
#' @param variables Number of marginals.
#' @param gfun Vector of inhomogeneity transforms.
#' @param gfun_pars List of parameters for the inhomogeneity functions.
#' @param scale Scale.
#'
#' @return An object of class \linkS4class{iph}.
#' @export
#'
#' @examples 
#' under_mph <- mph(structure = c("gcoxian", "general"), dimension = 4)
#' miph(under_mph, gfun = c("weibull", "pareto"), gfun_pars = list(c(2), c(3)))
miph <- function(mph = NULL, # object of class mPH
                 gfun = NULL, # vector of gfun for each marginal
                 gfun_pars = NULL, # List of gfun parameters for each marginal
                 alpha = NULL, # vector
                 S = NULL, # List
                 structure = NULL,
                 dimension = 3,
                 variables = NULL,
                 scale = 1) {
  if (all(is.null(c(gfun, gfun_pars)))) {
    stop("input inhomogeneity function and parameters")
  }
  d <- length(gfun)

  if (is.null(mph)) {
    mph <- mph(alpha = alpha, S = S, structure = structure, dimension = dimension, variables = d)
  }

  if (!all(gfun %in% c("pareto", "weibull", "lognormal", "loglogistic", "gompertz", "gev", "identity"))) {
    stop("invalid gfun for at least one marginal")
  }

  if (all(gfun %in% c("pareto", "weibull", "lognormal", "gompertz"))) {
    for (i in 1:d) {
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

  for (i in 1:d) {
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
  name <- if (is(mph, "miph")) mph@name else paste("inhomogeneous ", mph@name, sep = "")

  methods::new("miph",
    name = name,
    pars = mph@pars,
    gfun = list(
      name = gfun, # a vector
      pars = gfun_pars, # a list
      inverse = ginv, # a list
      inverse_prime = ginv_prime, # a list
      intensity = lambda, # a list
      intensity_prime = lambda_prime # a list
    ),
    scale = scale
  )
}

#' Show method for multivariate inhomogeneous phase-type distributions
#'
#' @param object An object of class \linkS4class{miph}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "miph", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@pars)
  cat("g-function name:", object@gfun$name, "\n")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@gfun$pars)
})

#' Simulation method for inhomogeneous multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{miph}.
#' @param n An integer of length of realization.
#'
#' @return A realization of independent and identically distributed inhomogeneous
#'  multivariate phase-type variables. If x is a MoE miph an array of dimension c(n,d,m) is returned,
#'  with d the number of marginals and m the number of initial distribution vectors.
#' @export
#'
#' @examples
#' under_mph <- mph(structure = c("general", "general"))
#' obj <- miph(under_mph, gfun = c("weibull", "pareto"), gfun_pars = list(c(2), c(3)))
#' sim(obj, 100)
setMethod("sim", c(x = "miph"), function(x, n = 1000) {
  name <- x@gfun$name
  pars <- x@gfun$pars
  scale <- x@scale

  if (is.vector(x@pars$alpha)) {
    U <- numeric(0)
    for (i in 1:length(name)) {
      if (name[i] %in% c("pareto", "weibull", "lognormal", "loglogistic", "gompertz")) {
        U <- cbind(U, scale * riph(n, name[i], x@pars$alpha, x@pars$S[[i]], pars[[i]]))
      }
      if (name[i] %in% c("gev")) {
        U <- cbind(U, scale * rmatrixgev(n, x@pars$alpha, x@pars$S[[i]], pars[[i]][1], pars[[i]][2], pars[[i]][3]))
      }
    }
  } else if (is.matrix(x@pars$alpha)) {
    U <- array(NA, dim = c(n, length(name), nrow(x@pars$alpha)))
    xx <- x
    for (m in 1:nrow(x@pars$alpha)) {
      U2 <- numeric(0)
      xx@pars$alpha <- x@pars$alpha[m, ]

      for (i in 1:length(name)) {
        if (name[i] %in% c("pareto", "weibull", "lognormal", "loglogistic", "gompertz")) {
          U2 <- cbind(U2, scale * riph(n, name[i], xx@pars$alpha, x@pars$S[[i]], pars[[i]]))
        }
        if (name[i] %in% c("gev")) {
          U2 <- cbind(U2, scale * rmatrixgev(n, xx@pars$alpha, x@pars$S[[i]], pars[[i]][1], pars[[i]][2], pars[[i]][3]))
        }
      }
      U[, , m] <- U2
    }
  }
  U
})

#' Marginal method for multivariate inhomogeneous phase-type distributions
#'
#' @param x An object of class \linkS4class{miph}.
#' @param mar Indicator of which marginal.
#' @return An object of the of class \linkS4class{iph}.
#' @export
#'
#' @examples
#' under_mph <- mph(structure = c("general", "general"))
#' obj <- miph(under_mph, gfun = c("weibull", "pareto"), gfun_pars = list(c(2), c(3)))
#' marginal(obj, 1)
setMethod("marginal", c(x = "miph"), function(x, mar = 1) {
  if (!(mar %in% 1:length(x@pars$S))) {
    stop("maringal provided not available")
  }
  S <- x@pars$S
  par <- x@gfun$pars
  iph(ph(alpha = x@pars$alpha, S = S[[mar]]), gfun = x@gfun$name[mar], gfun_pars = par[[mar]])
})

#' Density method for multivariate inhomogeneous phase-type distributions
#'
#' @param x An object of class \linkS4class{miph}.
#' @param delta Matrix with right-censoring indicators (1 uncensored, 0 right censored).
#' @param y A matrix of observations.
#'
#' @return A list containing the locations and corresponding density evaluations.
#' @export
#'
#' @examples
#' under_mph <- mph(structure = c("general", "general"))
#' obj <- miph(under_mph, gfun = c("weibull", "pareto"), gfun_pars = list(c(2), c(3)))
#' dens(obj, c(1, 2))
setMethod("dens", c(x = "miph"), function(x, y, delta = NULL) {
  alpha <- x@pars$alpha
  S <- x@pars$S
  gfun.pars <- x@gfun$pars
  d <- length(x@pars$S)

  if (is.matrix(y)) {
    n <- nrow(y)
  }
  if (is.vector(y)) {
    n <- 1
    y <- t(y)
  }

  if (length(delta) == 0) {
    delta <- matrix(1, nrow = n, ncol = d)
  }
  if (is.vector(delta)) {
    delta <- as.matrix(t(delta))
  }

  res <- numeric(n)

  if (is.vector(alpha)) { # if x regular miph
    p <- length(x@pars$alpha)

    for (j in 1:p) {
      in_vect <- rep(0, p)
      in_vect[j] <- 1
      aux <- matrix(NA, n, d)
      for (i in 1:d) {
        y_inv <- x@gfun$inverse[[i]](gfun.pars[[i]], y[, i])
        y_int <- x@gfun$intensity[[i]](gfun.pars[[i]], y[, i])
        for (m in 1:n) {
          if (delta[m, i] == 1) {
            aux[m, i] <- phdensity(y_inv[m], in_vect, S[[i]]) * y_int[m]
          } else {
            aux[m, i] <- phcdf(y_inv[m], in_vect, S[[i]], lower_tail = F)
          }
        }
      }
      res <- res + alpha[j] * apply(aux, 1, prod)
    }
  } else if (is.matrix(alpha)) { # if x is a mixture-of-expert miph
    p <- ncol(alpha)
    inter_res <- matrix(0, n, p)
    aux <- array(NA, c(n, d, p))

    for (j in 1:p) {
      in_vect <- rep(0, p)
      in_vect[j] <- 1

      for (i in 1:d) {
        y_inv <- x@gfun$inverse[[i]](gfun.pars[[i]], y[, i])
        y_int <- x@gfun$intensity[[i]](gfun.pars[[i]], y[, i])
        for (m in 1:n) {
          if (delta[m, i] == 1) {
            aux[m, i, j] <- phdensity(y_inv[m], in_vect, S[[i]]) * y_int[m]
          } else {
            aux[m, i, j] <- phcdf(y_inv[m], in_vect, S[[i]], lower_tail = F)
          }
        }
      }
    }

    for (m in 1:n) {
      for (j in 1:p) {
        inter_res[m, j] <- alpha[m, j] * prod(aux[m, , j])
      }
    }
    res <- rowSums(inter_res)
  }

 res
})

#' Distribution method for multivariate inhomogeneous phase-type distributions
#'
#' @param x An object of class \linkS4class{miph}.
#' @param y A matrix of observations.
#' @param lower.tail Logical parameter specifying whether lower tail (CDF) or
#'  upper tail is computed.
#'
#' @return A list containing the locations and corresponding CDF evaluations.
#' @export
#'
#' @examples
#' under_mph <- mph(structure = c("general", "general"))
#' obj <- miph(under_mph, gfun = c("weibull", "pareto"), gfun_pars = list(c(2), c(3)))
#' cdf(obj, c(1, 2))
setMethod("cdf", c(x = "miph"), function(x, y, lower.tail = TRUE) {
  alpha <- x@pars$alpha
  S <- x@pars$S
  gfun.pars <- x@gfun$pars
  d <- length(x@pars$S)

  if (is.matrix(y)) {
    n <- nrow(y)
  }
  if (is.vector(y)) {
    n <- 1
    y <- t(y)
  }

  res <- numeric(n)

  if (is.vector(alpha)) {
    p <- length(x@pars$alpha)

    for (j in 1:p) {
      in_vect <- rep(0, p)
      in_vect[j] <- 1
      aux <- matrix(NA, n, d)
      for (i in 1:d) {
        y_inv <- x@gfun$inverse[[i]](gfun.pars[[i]], y[, i])
        aux[, i] <- phcdf(y_inv, in_vect, S[[i]], lower.tail)
      }
      res <- res + alpha[j] * apply(aux, 1, prod)
    }
  } else if (is.matrix(alpha)) {
    p <- ncol(alpha)
    inter_res <- matrix(0, n, p)
    aux <- array(NA, c(n, d, p))

    for (j in 1:p) {
      in_vect <- rep(0, p)
      in_vect[j] <- 1

      for (i in 1:d) {
        y_inv <- x@gfun$inverse[[i]](gfun.pars[[i]], y[, i])
        y_int <- x@gfun$intensity[[i]](gfun.pars[[i]], y[, i])
        for (m in 1:n) {
          aux[m, i, j] <- phcdf(y_inv[m], in_vect, S[[i]], lower.tail)
        }
      }
    }

    for (m in 1:n) {
      for (j in 1:p) {
        inter_res[m, j] <- alpha[m, j] * prod(aux[m, , j])
      }
    }
    res <- rowSums(inter_res)
  }
  res
})
