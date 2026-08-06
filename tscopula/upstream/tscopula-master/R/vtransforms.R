#' Class of v-transforms
#'
#' This is the class of v-transforms. It contains the \linkS4class{VtransformI} subclass consisting of v-transforms
#' with an analytical expression for the inverse.
#'
#' @slot name a name for the v-transform of class character.
#' @slot Vtrans function to evaluate the v-transform.
#' @slot pars vector containing the named parameters of the v-transform.
#' @slot gradient function to evaluate the gradient of the v-transform.
#'
#' @export
#'
#' @examples
#' V2p(delta = 0.5, kappa = 1.2)
setClass("Vtransform", slots = list(
  name = "character", Vtrans = "function", pars = "numeric",
  gradient = "function"
))

#' Class of invertible v-transforms
#'
#' This class inherits from the \linkS4class{Vtransform} class and contains v-transforms
#' with an analytical expression for the inverse.
#'
#' @slot name a name for the v-transform of class character.
#' @slot Vtrans function to evaluate the v-transform.
#' @slot pars vector containing the named parameters of the v-transform.
#' @slot gradient function to evaluate the gradient of the v-transform.
#' @slot inverse function to evaluate the inverse of the v-transform.
#'
#' @export
#'
#' @examples
#' Vlinear(delta = 0.55)
setClass("VtransformI", contains = "Vtransform", slots = list(
  name = "character", Vtrans = "function",
  pars = "numeric", gradient = "function", inverse = "function"
))

#' Constructor function for symmetric v-transform
#'
#' @return An object of class \linkS4class{VtransformI}.
#' @export
#'
#' @examples
#' Vsymmetric()
Vsymmetric <- function() {
  new("VtransformI", name = "Vsymmetric", Vtrans = function(u) {
    abs(2 * u - 1)
  }, gradient = function(u) {
    slope <- rep(-2, length(u))
    slope[u > 0.5] <- 2
    slope
  }, inverse = function(v) {
    (1 - v) / 2
  })
}

#' Constructor function for degenerate v-transform
#'
#' @return An object of class \linkS4class{VtransformI}.
#' @export
#'
#' @examples
#' Vdegenerate()
Vdegenerate <- function() {
  new("VtransformI", name = "Vdegenerate", Vtrans = function(u) {
    u
  }, gradient = function(u) {
    rep(1, length(u))
  }, inverse = function(v) {
    v
  })
}

#' Constructor function for linear v-transform
#'
#' @param delta a value in (0, 1) specifying the fulcrum of the v-transform.
#'
#' @return An object of class \linkS4class{VtransformI}.
#' @export
#'
#' @examples
#' Vlinear(delta = 0.45)
Vlinear <- function(delta = 0.5) {
  new("VtransformI", name = "Vlinear", Vtrans = function(u, delta) {
    if (delta == 0) {
      return(u)
    } else {
      abs(u / delta - 1) * ((delta / (1 - delta))^(u > delta))
    }
  }, pars = c(delta = delta), gradient = function(u, delta) {
    slope <- rep(-1 / delta, length(u))
    slope[u > delta] <- 1 / (1 - delta)
    slope
  }, inverse = function(v, delta) {
    delta * (1 - v)
  })
}

#' Constructor function for 2-parameter v-transform
#'
#' @param delta a value in (0, 1) specifying the fulcrum of the v-transform.
#' @param kappa additional positive parameter of v-transform.
#'
#' @return An object of class \linkS4class{Vtransform}.
#' @export
#'
#' @examples
#' V2p(delta = 0.45, kappa = 1.2)
V2p <- function(delta = 0.5, kappa = 1) {
  new("Vtransform", name = "V2p", Vtrans = function(u, delta, kappa) {
    if (delta == 0)
      return(u)
    else if (delta == 1)
      return(1-u)
    else {
    ifelse(u <= delta, 1 - u - (1 - delta) * exp(-kappa * log(delta / u)), u - delta *
      exp(-(-log((1 - u) / (1 - delta)) / kappa)))
    }
  }, pars = c(delta = delta, kappa = kappa),
  gradient = function(u, delta, kappa) {
    slope <- rep(-1, length(u))
    arg1 <- log(delta / u[u <= delta])
    arg2 <- log((1 - delta) / (1 - u[u > delta]))
    slope[u <= delta] <- -1 - (1 - delta) * exp(-kappa * arg1) * kappa / u[u <= delta]
    slope[u > delta] <- 1 + delta * exp(-kappa * arg2) * kappa / (1 - u[u > delta])
    if (length(u[u == 0]) > 0) {
      if (kappa < 1) {
        val <- -Inf
      }
      if (kappa > 1) {
        val <- -1
      }
      if (kappa == 1) {
        val <- -1 / delta
      }
      slope[(u == 0)] <- val
    }
    if (length(u[u == 1]) > 0) {
      if (kappa > 1) {
        val2 <- Inf
      }
      if (kappa < 1) {
        val2 <- 1
      }
      if (kappa == 1) {
        val2 <- 1 / (1 - delta)
      }
      slope[(u == 1)] <- val2
    }
    slope
  })
}

#' Constructor function for 2-parameter beta v-transform
#'
#' @param delta a value in (0, 1) specifying the fulcrum of the v-transform.
#' @param kappa additional positive parameter of v-transform.
#'
#' @return An object of class \linkS4class{Vtransform}.
#' @export
#'
#' @examples
#' V2b(delta = 0.45, kappa = 1.2)
V2b <- function(delta = 0.5, kappa = 1) {
  new("Vtransform", name = "V2b", Vtrans = function(u, delta, kappa) {
    if (delta == 0)
      return(u)
    else if (delta == 1)
      return(1-u)
    else {
    suppressWarnings(ifelse(u <= delta, 1 - u - (1 - delta) * pbeta(
      u / delta, kappa,
      1 / kappa
    ), u - delta * qbeta((1 - u) / (1 - delta), kappa, 1 / kappa)))
    }
  }, pars = c(delta = delta, kappa = kappa),
  gradient = function(u, delta, kappa) {
    slope <- rep(-1, length(u))
    slope[u <= delta] <- -1 - (1 - delta) * dbeta(u[u <= delta] / delta, kappa, 1 / kappa) / delta
    slope[u > delta] <- 1 + delta / (dbeta(qbeta((1 - u[u > delta]) / (1 - delta), kappa, 1 / kappa), kappa, 1 / kappa) * (1 - delta))
    slope
  })
}

#' Constructor function for 3-parameter v-transform
#'
#' @param delta a value in (0, 1) specifying the fulcrum of the v-transform.
#' @param kappa additional positive parameter of v-transform.
#' @param xi additional positive parameter of v-transform.
#'
#' @return An object of class \linkS4class{Vtransform}.
#' @export
#'
#' @examples
#' V3p(delta = 0.45, kappa = 0.8, xi = 1.1)
V3p <- function(delta = 0.5, kappa = 1, xi = 1) {
  new("Vtransform", name = "V3p", Vtrans = function(u, delta, kappa, xi) {
    if (delta == 0)
      return(u)
    else if (delta == 1)
      return(1-u)
    else {
    ifelse(u <= delta, 1 - u - (1 - delta) * exp(-kappa * (log(delta / u))^xi), u -
      delta * exp(-(-log((1 - u) / (1 - delta)) / kappa)^(1 / xi)))
    }
  }, pars = c(delta = delta, kappa = kappa, xi = xi),
  gradient = function(u, delta, kappa, xi) {
    slope <- rep(-1, length(u))
    arg1 <- log(delta / u[u <= delta])
    arg2 <- log((1 - delta) / (1 - u[u > delta]))
    slope[u <= delta] <- -1 - (1 - delta) * exp(-kappa * arg1^xi) * arg1^(xi -
      1) * (xi * kappa / u[u <= delta])
    slope[u > delta] <- 1 + delta * exp(-(kappa * arg2)^(1 / xi)) * arg2^(1 / xi - 1) *
      kappa^(1 / xi) / (xi * (1 - u[u > delta]))
    if (length(u[u == 0]) > 0) {
      if ((xi < 1) || ((xi == 1) && (kappa < 1))) {
        val <- -Inf
      }
      if ((xi > 1) || ((xi == 1) && (kappa > 1))) {
        val <- -1
      }
      if ((xi == 1) && (kappa == 1)) {
        val <- -1 / delta
      }
      slope[(u == 0)] <- val
    }
    if (length(u[u == 1]) > 0) {
      if ((xi > 1) || ((xi == 1) && (kappa < 1))) {
        val2 <- Inf
      }
      if ((xi < 1) || ((xi == 1) && (kappa > 1))) {
        val2 <- 1
      }
      if ((xi == 1) && (kappa == 1)) {
        val2 <- 1 / (1 - delta)
      }
      slope[(u == 1)] <- val2
    }
    slope
  })
}

#' Constructor function for 3-parameter beta v-transform
#'
#' @param delta a value in (0, 1) specifying the fulcrum of the v-transform.
#' @param kappa additional positive parameter of v-transform.
#' @param xi additional positive parameter of v-transform.
#'
#' @return An object of class \linkS4class{Vtransform}.
#' @export
#'
#' @examples
#' V3b(delta = 0.45, kappa = 1.2, xi = 1.2)
V3b <- function(delta = 0.5, kappa = 1, xi = 1) {
  new("Vtransform", name = "V2b", Vtrans = function(u, delta, kappa, xi) {
    if (delta == 0)
      return(u)
    else if (delta == 1)
      return(1-u)
    else {
    suppressWarnings(ifelse(u <= delta, 1 - u - (1 - delta) * pbeta(
      u / delta, kappa,
      xi
    ), u - delta * qbeta((1 - u) / (1 - delta), kappa, xi)))
    }
  }, pars = c(delta = delta, kappa = kappa, xi = xi),
  gradient = function(u, delta, kappa, xi) {
    slope <- rep(-1, length(u))
    slope[u <= delta] <- -1 - (1 - delta) * dbeta(u[u <= delta] / delta, kappa, xi) / delta
    slope[u > delta] <- 1 + delta / (dbeta(qbeta((1 - u[u > delta]) / (1 - delta), kappa, xi), kappa, xi) * (1 - delta))
    slope
  })
}

#' Evaluate a v-transform
#'
#' @param x an object of class \linkS4class{Vtransform}.
#' @param u a vector or time series with values in [0, 1].
#'
#'
#' @return A vector or time series with values in [0, 1].
#' @export
#'
#' @examples
#' vtrans(Vsymmetric(), c(0, 0.25, 0.5, 0.75, 1))
vtrans <- function(x, u) {
  do.call(x@Vtrans, append(x@pars, list(u = u)))
}

#' Calculate gradient of v-transform
#'
#' @param x an object of class \linkS4class{Vtransform}.
#' @param u a vector or time series with values in [0, 1].
#'
#' @return A vector or time series of values of gradient.
#' @export
#'
#' @examples
#' vgradient(Vsymmetric(), c(0, 0.25, 0.5, 0.75, 1))
vgradient <- function(x, u) {
  do.call(x@gradient, append(x@pars, list(u = u)))
}

#' Calculate inverse of v-transform
#'
#' If the \linkS4class{Vtransform} object is also a \linkS4class{VtransformI} object (an
#' invertible v-transform) then the analytical inverse is used. Otherwise
#' an inverse is found by numerical root finding with \code{\link[stats]{uniroot}}.
#'
#' @param x an object ofc lass \linkS4class{Vtransform}.
#' @param v a vector or time series with values in [0, 1].
#' @param tol the desired accuracy (convergence tolerance) that is passed to
#' \code{uniroot} if numerical inversion is used.
#'
#' @return A vector or time series with values in [0, 1].
#' @export
#'
#' @examples
#' vinverse(Vsymmetric(), c(0, 0.25, 0.5, 0.75, 1))
vinverse <- function(x, v, tol = .Machine$double.eps^0.75) {
  if (is(x,"VtransformI")) {
    do.call(x@inverse, append(x@pars, list(v = v)))
  } else {
    vecinverse <- Vectorize(function(v, vfunc, pars, tol) {
      uniroot(function(t) {
        do.call(vfunc, append(pars, list(u = t))) - v
      }, lower = 0, upper = pars["delta"], tol = tol)$root
    }, "v")
    vecinverse(v, x@Vtrans, x@pars, tol)
  }
}

#' Calculate conditional down probability of v-transform
#'
#' @param x an object of class \linkS4class{Vtransform}.
#' @param v a vector or time series with values in [0, 1].
#'
#' @return A vector or time series of values of gradient.
#' @export
#'
#' @examples
#' vdownprob(V2p(delta = 0.55, kapp = 1.2), c(0, 0.25, 0.5, 0.75, 1))
vdownprob <- function(x, v) {
  -1 / vgradient(x, vinverse(x, v))
}

#' Stochastic inverse of a v-transform
#'
#' @param x an object of class \linkS4class{Vtransform}.
#' @param v a vector, matrix or time series with values in [0, 1].
#' @param tscopula a time series copula object.
#' @param tol the desired accuracy (convergence tolerance) that is passed to
#' \code{uniroot} if numerical inversion is used.
#'
#' @return A vector, matrix or time series with values in [0, 1].
#' @export
#'
#'
#' @examples
#' stochinverse(Vsymmetric(), c(0, 0.25, 0.5, 0.75, 1))
stochinverse <- function(x, v, tscopula = NULL, tol = .Machine$double.eps^0.75) {
  vinv <- vinverse(x, v, tol)
  pdown <- -1 / vgradient(x, vinv)
  if (!(is.null(tscopula))) {
    W <- sim(tscopula, length(v))
  } else {
    W <- runif(length(v))
  }
  output <- ifelse(W <= pdown, vinv, v + vinv)
  if (!(is.null(attributes(v)))) {
    attributes(output) <- attributes(v)
  }
  output
}

#' Plot method for Vtransform class
#'
#' Plots the v-transform as well as its gradient or inverse. Can also plot the
#' conditional probability that a series PIT falls below the fulcrum for a
#' given volatility PIT value v.
#'
#' @param x an object of class \linkS4class{Vtransform}.
#' @param type type of plot: 'transform' for plot of transform, 'inverse' for plot of inverse,
#' 'gradient' for plot of gradient or 'pdown' for plot of conditional probability.
#' @param shading logical variable specifying whether inadmissible zone for v-transform
#' should be shaded
#' @param npoints number of plotting points along x-axis.
#' @param lower the lower x-axis value for plotting.
#' @param upper the upper x-axis value for plotting
#'
#' @return No return value, generates plot.
#' @export
#'
#'
#' @examples
#' plot(Vsymmetric())
#' plot(V2p(delta = 0.45, kappa = 0.8), type = "inverse")
#' plot(V2p(delta = 0.45, kappa = 0.8), type = "gradient")
setMethod("plot", c(x = "Vtransform", y = "missing"), function(x, type = "transform",
                                                               shading = TRUE, npoints = 200, lower = 0, upper = 1) {
  delta <- ifelse(is.element("delta", names(x@pars)), x@pars["delta"], 0.5)
  switch(type, inverse = {
    vvals <- seq(from = max(lower, 0), to = min(upper, 1), length = npoints)
    plot(vvals, vinverse(x, vvals), xlab = "v", ylab = "Vinv(v)", type = "l")
  }, gradient = {
    uvals <- seq(from = max(lower, 0), to = min(upper, 1), length = npoints)
    plot(uvals, vgradient(x, uvals), xlab = "u", ylab = "Vprime(u)", type = "l")
  }, pdown = {
    vvals <- seq(from = max(lower, 0), to = min(upper, 1), length = npoints)
    plot(vvals, vdownprob(x, vvals), xlab = "v", ylab = "Delta(v)", type = "l")
  }, transform = {
    uvals <- seq(from = max(lower, 0), to = min(upper, 1), length = npoints)
    if ((delta > lower) & (delta < upper)) uvals <- sort(c(uvals, delta))
    plot(uvals, vtrans(x, uvals), xlab = "u", ylab = "V(u)", type = "l")
    if (shading) {
      # colchoice = 'gray97'
      colchoice <- "gray90"
      polygon(c(0, 0, delta), c(delta, 0, 0), col = colchoice, border = NA)
      polygon(c(delta, 1, 1), c(0, 0, 1 - delta), col = colchoice, border = NA)
      polygon(c(0, delta, delta), c(1, 1, 1 - delta), col = colchoice, border = NA)
      polygon(c(delta, delta, 1), c(delta, 1, 1), col = colchoice, border = NA)
    }
  }, stop("Not a plot method for v-transform."))
})

#' Compute coincidence probability for v-transform
#'
#' Computes the probability that if we v-transform a uniform
#' random variable and then stochastically invert the
#' v-transform, we get back to the original value.
#'
#' @param x an object of class \linkS4class{Vtransform}.
#'
#' @return The probability of coincidence.
#' @export
#'
#' @examples
#' pcoincide(Vlinear(delta = 0.4))
#' pcoincide(V3p(delta = 0.45, kappa = 0.5, xi = 1.3))
pcoincide <- function(x) {
  if (x@name == "symmetric") {
    return(0.5)
  }
  integrand <- function(v, Vtransform) (vdownprob(Vtransform, v) - Vtransform@pars[1])^2
  delta <- x@pars[1]
  varDelta <- integrate(integrand, 0, 1, Vtransform = x)$value
  unname(delta^2 + (1 - delta)^2 + 2 * varDelta)
}

#' @describeIn Vtransform Show method for Vtransform class
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("show", "Vtransform", function(object) {
  cat("name: ", object@name, "\n", "parameters: ",
    "\n",
    sep = ""
  )
  print(object@pars)
})

#' @describeIn Vtransform Coef method for Vtransform class
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("coef", "Vtransform", function(object) {
  object@pars
})
