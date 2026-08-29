#' Bivariate phase-type distributions
#'
#' Class of objects for bivariate phase-type distributions.
#'
#' @slot name Name of the phase-type distribution.
#' @slot pars A list comprising of the parameters.
#' @slot fit A list containing estimation information.
#'
#' @return Class object.
#' @export
#'
setClass("bivph",
  slots = list(
    name = "character",
    pars = "list",
    fit = "list"
  ),
  prototype = list(
    name = NA_character_,
    pars = list(),
    fit = list()
  )
)

#' Constructor function for bivariate phase-type distributions
#'
#' @param alpha A probability vector.
#' @param S11 A sub-intensity matrix.
#' @param S12 A matrix.
#' @param S22 A sub-intensity matrix.
#' @param dimensions The dimensions of the bivariate phase-type (if no parameters are provided).
#'
#' @return An object of class \linkS4class{bivph}.
#' @export
#'
#' @examples
#' bivph(dimensions = c(3, 3))
#' S11 <- matrix(c(-1, .5, .5, -1), 2, 2)
#' S12 <- matrix(c(.2, .4, .3, .1), 2, 2)
#' S22 <- matrix(c(-2, 0, 1, -1), 2, 2)
#' bivph(alpha = c(.5, .5), S11, S12, S22)
bivph <- function(alpha = NULL, S11 = NULL, S12 = NULL, S22 = NULL, dimensions = c(3, 3)) {
  if (any(is.null(alpha)) & any(is.null(S11))) {
    rs <- random_structure_bivph(dimensions[1], dimensions[2])
    alpha <- rs[[1]]
    S11 <- rs[[2]]
    S12 <- rs[[3]]
    S22 <- rs[[4]]
    name <- "random"
  } else {
    if (dim(S11)[1] != dim(S11)[2]) {
      stop("matrix S11 should be square")
    } else if (dim(S22)[1] != dim(S22)[2]) {
      stop("matrix S22 should be square")
    } else if (dim(S11)[1] != dim(S12)[1]) {
      stop("incompatible dimensions of S11 and S12")
    } else if (dim(S22)[1] != dim(S12)[2]) {
      stop("incompatible dimensions of S12 and S22")
    } else if (length(alpha) != dim(S11)[1]) {
      stop("incompatible dimensions of alpha and S11")
    }
    name <- "custom"
  }
  methods::new("bivph",
    name = paste(name, " bivph(", length(alpha), ")", sep = ""),
    pars = list(alpha = alpha, S11 = S11, S12 = S12, S22 = S22)
  )
}

#' Show method for bivariate phase-type distributions
#'
#' @param object An object of class \linkS4class{bivph}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "bivph", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@pars)
})

#' Density method for bivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{bivph}.
#' @param y A matrix of locations.
#'
#' @return A vector containing the joint density evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' dens(obj, matrix(c(0.5, 1), ncol = 2))
setMethod("dens", c(x = "bivph"), function(x, y) {
  if (is.vector(y)) {
    y <- t(y)
  }
  bivph_density(y, x@pars$alpha, x@pars$S11, x@pars$S12, x@pars$S22)
})

#' Simulation method for bivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{bivph}.
#' @param n An integer of length of realization.
#'
#' @return A realization of independent and identically distributed bivariate
#'  phase-type vector.
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' sim(obj, n = 100)
setMethod("sim", c(x = "bivph"), function(x, n = 1000) {
  p1_aux <- dim(x@pars$S11)[1]
  p2_aux <- dim(x@pars$S22)[1]
  alpha_aux <- c(x@pars$alpha, rep(0, p2_aux))
  S_aux <- merge_matrices(x@pars$S11, x@pars$S12, x@pars$S22)
  R_aux <- matrix(c(c(rep(1, p1_aux), rep(0, p2_aux)), c(rep(0, p1_aux), rep(1, p2_aux))), ncol = 2)
  rMPHstar(n, alpha_aux, S_aux, R_aux)
})

#' Coef method for bivph class
#'
#' @param object An object of class \linkS4class{bivph}.
#'
#' @return Parameters of bivariate phase-type model.
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' coef(obj)
setMethod("coef", c(object = "bivph"), function(object) {
  object@pars
})

#' Moment method for bivph class
#'
#' @param x An object of class \linkS4class{bivph}.
#' @param k A vector with the location.
#' @return An real value.
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' moment(obj, c(1, 1))
setMethod("moment", c(x = "bivph"), function(x, k = c(1, 1)) {
  if (all(k == 0) | any(k < 0)) {
    stop("k should be non-negative and not zero")
  }
  if (any((k %% 1) != 0)) {
    stop("k should be an integer")
  }
  if (methods::is(x, "biviph")) {
    warning("moment of undelying bivph structure is provided for biviph objects")
  }
  ee <- rep(1, nrow(x@pars$S22))
  return(factorial(k[1]) * factorial(k[2]) * x@pars$alpha %*% matrix_power(k[1] + 1, base::solve(-x@pars$S11)) %*% x@pars$S12 %*% matrix_power(k[2], base::solve(-x@pars$S22)) %*% ee)
})

#' Mean Method for bivph class
#'
#' @param x An object of class \linkS4class{bivph}.
#'
#' @return The mean of the bivariate phase-type distribution.
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' mean(obj)
setMethod("mean", c(x = "bivph"), function(x) {
  if (methods::is(x, "biviph")) {
    warning("mean of undelying bivph structure is provided for biviph objects")
  }
  suppressWarnings(c(moment(x, c(1, 0)), moment(x, c(0, 1))))
})

#' Var method for bivph class
#'
#' @param x An object of class \linkS4class{bivph}.
#'
#' @return The covariance matrix of the bivariate phase-type distribution.
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' var(obj)
setMethod("var", c(x = "bivph"), function(x) {
  if (methods::is(x, "biviph")) {
    warning("covariance matrix of undelying bivph structure is provided for biviph objects")
  }
  re <- matrix(0, 2, 2)
  re[1, 1] <- suppressWarnings(moment(x, c(2, 0)) - moment(x, c(1, 0))^2)
  re[1, 2] <- suppressWarnings(moment(x, c(1, 1)) - moment(x, c(1, 0)) * moment(x, c(0, 1)))
  re[2, 1] <- re[1, 2]
  re[2, 2] <- suppressWarnings(moment(x, c(0, 2)) - moment(x, c(0, 1))^2)
  re
})

#' Cor method for bivph class
#'
#' @param x An object of class \linkS4class{bivph}.
#'
#' @return The correlation matrix of the bivariate phase-type distribution.
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' cor(obj)
setMethod("cor", c(x = "bivph"), function(x) {
  if (methods::is(x, "biviph")) {
    warning("correlation matrix of undelying bivph structure is provided for biviph objects")
  }
  suppressWarnings(stats::cov2cor(var(x)))
})

#' Laplace method for bivph class
#'
#' @param x An object of class \linkS4class{mph}.
#' @param r A matrix of real values.
#'
#' @return A vector containing the corresponding Laplace transform evaluations.
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' laplace(obj, matrix(c(0.5, 1), ncol = 2))
setMethod("laplace", c(x = "bivph"), function(x, r) {
  if (methods::is(x, "biviph")) {
    warning("Laplace transform of undelying bivph structure is provided for biviph objects")
  }

  if (is.vector(r)) {
    n <- 1
    r <- t(r)
  }

  lim1 <- max(Re(eigen(x@pars$S11)$values))
  lim2 <- max(Re(eigen(x@pars$S22)$values))
  if (any(r[, 1] <= lim1) | any(r[, 2] <= lim2)) {
    stop("r should be above the largest real eigenvalue of S")
  }

  bivph_laplace(r, x@pars$alpha, x@pars$S11, x@pars$S12, x@pars$S22)
})

#' Mgf method for bivph class
#'
#' @param x An object of class \linkS4class{mph}.
#' @param r A matrix of real values.
#'
#' @return A vector containing the corresponding mgf evaluations.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- bivph(dimensions = c(3, 3))
#' mgf(obj, matrix(c(0.5, 0.1), ncol = 2))
setMethod("mgf", c(x = "bivph"), function(x, r) {
  if (methods::is(x, "biviph")) {
    warning("mgf of undelying bivph structure is provided for biviph objects")
  }

  if (is.vector(r)) {
    n <- 1
    r <- t(r)
  }

  lim1 <- -max(Re(eigen(x@pars$S11)$values))
  lim2 <- -max(Re(eigen(x@pars$S22)$values))
  if (any(r[, 1] > lim1) | any(r[, 2] > lim2)) {
    stop("r should be below the negative largest real eigenvalue of S")
  }

  bivph_laplace(-r, x@pars$alpha, x@pars$S11, x@pars$S12, x@pars$S22)
})

#' Marginal method for bivph class
#'
#' @param x An object of class \linkS4class{bivph}.
#' @param mar Indicator of which marginal.
#' @return An object of the of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' marginal(obj, 1)
setMethod("marginal", c(x = "bivph"), function(x, mar = 1) {
  if (!(mar %in% 1:2)) {
    stop("maringal provided not available")
  }
  if (mar == 1) {
    x0 <- ph(alpha = x@pars$alpha, S = x@pars$S11)
  } else {
    alpha0 <- x@pars$alpha %*% base::solve(-x@pars$S11) %*% x@pars$S12
    x0 <- ph(alpha = alpha0, S = x@pars$S22)
  }
  x0
})

#' Linear combination method for bivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{bivph}.
#' @param w A vector with non-negative entries.
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' linCom(obj, c(1, 0))
setMethod("linCom", c(x = "bivph"), function(x, w = c(1, 1)) {
  if (length(w) != 2) {
    stop("vector of wrong dimension")
  }
  if (any(w < 0)) {
    stop("vector with negative entries")
  }
  if (all(w == 0)) {
    stop("vector with all entries zero")
  }
  if (methods::is(x, "biviph")) {
    warning("Linear combination of undelying bivph structure is provided for biviph objects")
  }
  p1_aux <- dim(x@pars$S11)[1]
  p2_aux <- dim(x@pars$S22)[1]
  alpha_aux <- c(x@pars$alpha, rep(0, p2_aux))
  S_aux <- merge_matrices(x@pars$S11, x@pars$S12, x@pars$S22)
  R_aux <- matrix(c(c(rep(1, p1_aux), rep(0, p2_aux)), c(rep(0, p1_aux), rep(1, p2_aux))), ncol = 2)
  L <- linear_combination(w, alpha_aux, S_aux, R_aux)
  ph(alpha = L$alpha, S = L$S)
})

#' Fit method for bivph Class
#'
#' @param x An object of class \linkS4class{bivph}.
#' @param y A matrix with the data.
#' @param weight Vector of weights.
#' @param stepsEM Number of EM steps to be performed.
#' @param every Number of iterations between likelihood display updates.
#' @param maxit Maximum number of iterations when optimizing g functions.
#' @param reltol Relative tolerance when optimizing g functions.
#'
#' @return An object of class \linkS4class{bivph}.
#'
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' data <- sim(obj, n = 100)
#' fit(obj, data, stepsEM = 100, every = 50)
setMethod(
  "fit", c(x = "bivph"),
  function(x,
           y,
           weight = numeric(0),
           stepsEM = 1000,
           maxit = 100,
           reltol = 1e-8,
           every = 10) {
    if (!all(y > 0)) {
      stop("data should be positive")
    }
    if (!all(weight >= 0)) {
      stop("weights should be non-negative")
    }
    if (length(weight) == 0) {
      weight <- rep(1, length(y[, 1]))
    }

    is_biviph <- methods::is(x, "biviph")

    bivph_par <- x@pars
    alpha_fit <- clone_vector(bivph_par$alpha)
    S11_fit <- clone_matrix(bivph_par$S11)
    S12_fit <- clone_matrix(bivph_par$S12)
    S22_fit <- clone_matrix(bivph_par$S22)

    options(digits.secs = 4)
    cat(format(Sys.time(), format = "%H:%M:%OS"), ": EM started", sep = "")
    cat("\n", sep = "")

    if (!is_biviph) {
      for (k in 1:stepsEM) {
        EMstep_bivph(alpha_fit, S11_fit, S12_fit, S22_fit, y, weight)
        if (k %% every == 0) {
          cat("\r", "iteration:", k,
            ", logLik:", logLikelihoodbivPH(alpha_fit, S11_fit, S12_fit, S22_fit, y, weight),
            sep = " "
          )
        }
      }
      x@pars$alpha <- alpha_fit
      x@pars$S11 <- S11_fit
      x@pars$S12 <- S12_fit
      x@pars$S22 <- S22_fit
    } else if (is_biviph) {
      par_name <- x@gfun$name
      par_g <- x@gfun$pars
      inv_g <- x@gfun$inverse
      opt_fun <- biviph_LL

      for (k in 1:stepsEM) {
        trans <- clone_matrix(y)
        for (i in 1:2) {
          if (x@gfun$name[i] != "gev") {
            trans[, i] <- inv_g[[i]](par_g[[i]], y[, i])
          } else {
            t <- inv_g[[i]](par_g[[i]], y[, i], rep(1, nrow(y)))
            trans[, i] <- t$obs
          }
        }

        EMstep_bivph(alpha_fit, S11_fit, S12_fit, S22_fit, trans, weight)

        x@pars$alpha <- alpha_fit
        x@pars$S11 <- S11_fit
        x@pars$S12 <- S12_fit
        x@pars$S22 <- S22_fit

        opt <- suppressWarnings(
          stats::optim(
            par = par_g,
            fn = opt_fun,
            x = x,
            obs = y,
            hessian = F,
            control = list(
              maxit = maxit,
              reltol = reltol,
              fnscale = -1
            )
          )
        )

        par_g <- as.list(opt$par)
        if (k %% every == 0) {
          cat("\r", ", iteration:", k,
            ", logLik:", opt$value,
            sep = " "
          )
        }
      }
      x@gfun$pars <- par_g
    }

    cat("\n", format(Sys.time(), format = "%H:%M:%OS"), ": EM finalized", sep = "")
    cat("\n", sep = "")

    x
  }
)

# multivariate loglikelihood to be optimized
biviph_LL <- function(x,
                      obs,
                      gfun_pars) {
  x@gfun$pars <- gfun_pars
  res <- dens(x = x, y = obs)
  sum(log(res))
}
