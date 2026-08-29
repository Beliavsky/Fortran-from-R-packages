#' Survival analysis for phase-type distributions
#'
#' Class of objects for inhomogeneous phase-type distributions
#'
#' @slot coefs Coefficients of the survival regression object.
#' @slot type Type of survival object.
#'
#' @return Class object
#' @export
#'
setClass("sph",
  contains = c("iph"),
  slots = list(
    coefs = "list",
    type = "character"
  )
)

#' Constructor function for survival phase-type objects
#'
#' @param x An object of class \linkS4class{ph}.
#' @param coefs Coefficients of the survival regression object.
#' @param type Type of survival object.
#'
#' @return An object of class \linkS4class{sph}.
#' @export
#'
sph <- function(x = NULL, coefs = list(B = numeric(0), C = numeric(0)), type = "reg") {
  if (!type %in% c("reg", "reg2", "aft")) {
    stop("type must be one of : reg, reg2, aft")
  }
  if (!is(x, "iph")) {
    if (!is(x, "ph")) {
      ph <- ph(structure = "general")
    }
    gfun <- list(name = "identity", pars = numeric(0))
  } else {
    gfun <- x@gfun
  }
  new("sph",
    name = paste("survival", type, x@name, sep = " "),
    pars = x@pars,
    gfun = gfun,
    coefs = coefs,
    type = type,
    scale = 1,
    fit = x@fit
  )
}

#' Show method for survival phase-type objects
#'
#' @param object An object of class \linkS4class{sph}.
#' @export
#'
setMethod("show", "sph", function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@pars)
  cat("g-function name: ", object@gfun$name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@gfun$pars)
  cat("coefficients: ", "\n", sep = "")
  print(object@coefs)
})

#' Coef method for sph Class
#'
#' @param object An object of class \linkS4class{sph}.
#'
#' @return Parameters of sph model.
#' @export
#'
setMethod("coef", c(object = "sph"), function(object) {
  L <- append(object@pars, unname(object@gfun$pars))
  names(L)[3:(2 + length(object@gfun$pars))] <- names(object@gfun$pars)
  append(L, object@coefs)
})

#' Evaluation method for sph Class
#'
#' @param x An object of class \linkS4class{sph}.
#' @param subject Covariates of a single subject.
#'
#' @return A \linkS4class{ph} model.
#' @export
#'
setMethod("evaluate", c(x = "sph"), function(x, subject) {
  if (x@gfun$name == "identity") {
    if (x@type == "aft") {
      z <- ph(
        alpha = x@pars$alpha,
        S = x@pars$S / as.numeric(exp(subject %*% x@coefs$B))
      )
    } else if (x@type == "reg") {
      z <- ph(
        alpha = x@pars$alpha,
        S = as.numeric(exp(subject %*% x@coefs$B)) * x@pars$S
      )
    }
  } else {
    if (x@type == "aft") {
      z <- iph(
        gfun = x@gfun$name,
        gfun_pars = x@gfun$pars,
        alpha = x@pars$alpha,
        S = x@pars$S,
        scale = as.numeric(exp(subject %*% x@coefs$B))
      )
    } else if (x@type == "reg") {
      z <- iph(
        gfun = x@gfun$name,
        gfun_pars = x@gfun$pars,
        alpha = x@pars$alpha,
        S = as.numeric(exp(subject %*% x@coefs$B)) * x@pars$S,
      )
    } else if (x@type == "reg2") {
      z <- iph(
        gfun = x@gfun$name,
        gfun_pars = x@gfun$pars * as.numeric(exp(subject %*% x@coefs$C)),
        alpha = x@pars$alpha,
        S = as.numeric(exp(subject %*% x@coefs$B)) * x@pars$S,
      )
    }
  }
  z
})

#' Fisher information method for sph class
#'
#' @param x An object of class \linkS4class{sph}.
#' @param y Independent variate.
#' @param X Matrix of covariates.
#' @param w Weights.
#'
#' @return A matrix.
#' @export
#'
setMethod("Fisher", c(x = "sph"), function(x, y, X, w = numeric(0)) {
  if (x@type != "reg") {
    stop("method not implemented")
  }
  if (length(w) == 0) w <- rep(1, length(y))
  n <- ncol(X) + 1
  result <- matrix(0, n, n)
  for (j in 1:length(w)) {
    v <- unlist(log_lik_derivative(x, z = y[j], Z = X[j, ]))
    result <- result + w[j] * outer(v, v)
  }
  result
})

log_lik_derivative <- function(x, z, Z) {
  p <- length(Z) + 1
  result <- numeric(p)

  gfn <- x@gfun$inverse
  gfn_pr <- x@gfun$pars
  gf <- gfn(gfn_pr, z)
  alpha <- x@pars$alpha
  ee <- rep(1, length(alpha))
  S <- x@pars$S
  B <- x@coefs$B
  M <- matrix_exponential(exp(sum(Z * B)) * gf * S)
  ratio <- x@gfun$intensity_prime(gfn_pr, z) / x@gfun$intensity(gfn_pr, z)
  a <- -alpha %*% M %*% S %*% S %*% ee
  b <- -alpha %*% M %*% S %*% ee
  result[1] <- ratio + a * exp(sum(Z * B)) * x@gfun$inverse_prime(gfn_pr, z) / b
  for (i in 2:p) {
    a2 <- Z[i - 1] * gf * exp(sum(Z * B)) * a
    result[i] <- Z[i - 1] + a2 / b
  }
  result
}
