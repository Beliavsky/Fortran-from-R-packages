#' Multivariate phase-type distributions
#'
#' Class of objects for multivariate phase-type distributions.
#'
#' @slot name Name of the phase type distribution.
#' @slot pars A list comprising of the parameters.
#' @slot fit A list containing estimation information.
#'
#' @return Class object.
#' @export
#'
setClass("mph",
  slots = list(
    name = "character",
    pars = "list",
    fit = "list"
  )
)

#' Constructor function for multivariate phase-type distributions
#'
#' @param alpha A probability vector.
#' @param S A list of sub-intensity matrices.
#' @param structure A vector of valid ph structures.
#' @param dimension The dimension of the ph structure (if provided).
#' @param variables The dimension of the multivariate phase-type.
#'
#' @return An object of class \linkS4class{mph}.
#' @export
#'
#' @examples
#' mph(structure = c("gcoxian", "general"), dimension = 5)
mph <- function(alpha = NULL, S = NULL, structure = NULL, dimension = 3, variables = NULL) {
  if (any(is.null(alpha)) & any(is.null(S)) & is.null(structure)) {
    stop("input a vector and matrix, or a structure")
  }
  if (is.null(variables)) {
    variables <- length(structure)
  }
  if (!any(is.null(structure))) {
    rs <- random_structure(dimension, structure = structure[1])
    alpha <- rs[[1]]
    S <- list()
    S[[1]] <- rs[[2]]
    for (i in 2:variables) {
      S[[i]] <- random_structure(dimension, structure = structure[i])[[2]]
    }
    name <- structure
  } else {
    name <- "custom"
  }
  methods::new("mph",
    name = paste(name, " mph(", length(alpha), ")", sep = " "),
    pars = list(alpha = alpha, S = S)
  )
}

#' Show method for multivariate phase-type distributions
#'
#' @param object An object of class \linkS4class{mph}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "mph", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@pars)
  cat("number of variables: ", length(object@pars$S), "\n", sep = "")
})

#' Simulation method for multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{mph}.
#' @param n Length of realization.
#' @param equal_marginals Non-negative integer. If positive, it specifies
#' the number of marginals to simulate from, all from the first matrix.
#'
#' @return A realization of a multivariate phase-type distribution.
#' @export
#'
#' @examples
#' obj <- mph(structure = c("general", "general"))
#' sim(obj, 100)
setMethod("sim", c(x = "mph"), function(x, n = 1000, equal_marginals = 0) {
  if (is.vector(x@pars$alpha)) p <- length(x@pars$alpha)
  if (is.matrix(x@pars$alpha)) p <- ncol(x@pars$alpha)

  if (equal_marginals == 0) {
    d <- length(x@pars$S)
    states <- 1:p
    result <- matrix(NA, n, d)
    for (i in 1:n) {
      state <- sample(states, 1, prob = x@pars$alpha)
      in_vect <- rep(0, p)
      in_vect[state] <- 1
      for (j in 1:d) {
        result[i, j] <- rphasetype(1, in_vect, x@pars$S[[j]])
      }
    }
  } else {
    d <- equal_marginals
    states <- 1:p
    result <- matrix(NA, n, d)
    for (i in 1:n) {
      state <- sample(states, 1, prob = x@pars$alpha)
      in_vect <- rep(0, p)
      in_vect[state] <- 1
      for (j in 1:d) {
        result[i, j] <- rphasetype(1, in_vect, x@pars$S[[1]])
      }
    }
  }
  result
})

#' Marginal method for multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{mph}.
#' @param mar Indicator of which marginal.
#' @return An object of the of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' obj <- mph(structure = c("general", "general"))
#' marginal(obj, 1)
setMethod("marginal", c(x = "mph"), function(x, mar = 1) {
  if (!(mar %in% 1:length(x@pars$S))) {
    stop("maringal provided not available")
  }
  S <- x@pars$S
  ph(alpha = x@pars$alpha, S = S[[mar]])
})

#' Density method for multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{mph}.
#' @param delta Matrix with right-censoring indicators (1 uncensored, 0 right censored).
#' @param y A matrix of observations.
#'
#' @return A list containing the locations and corresponding density evaluations.
#' @export
#'
#' @examples
#' obj <- mph(structure = c("general", "general"))
#' dens(obj, matrix(c(0.5, 1), ncol = 2))
setMethod("dens", c(x = "mph"), function(x, y, delta = NULL) {
  alpha <- x@pars$alpha
  S <- x@pars$S
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

  if (is.vector(alpha)) {
    p <- length(x@pars$alpha)

    for (j in 1:p) {
      in_vect <- rep(0, p)
      in_vect[j] <- 1
      aux <- matrix(NA, n, d)
      for (i in 1:d) {
        for (m in 1:n) {
          if (delta[m, i] == 1) {
            aux[m, i] <- phdensity(y[m, i], in_vect, S[[i]])
          } else {
            aux[m, i] <- phcdf(y[m, i], in_vect, S[[i]], lower_tail = F)
          }
        }
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
        for (m in 1:n) {
          if (delta[m, i] == 1) {
            aux[m, i, j] <- phdensity(y[m, i], in_vect, S[[i]])
          } else {
            aux[m, i, j] <- phcdf(y[m, i], in_vect, S[[i]], lower_tail = F)
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

#' Distribution method for multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{mph}.
#' @param y A matrix of observations.
#' @param lower.tail Logical parameter specifying whether lower tail (CDF) or
#'  upper tail is computed.
#'
#' @return A list containing the locations and corresponding CDF evaluations.
#' @export
#'
#' @examples
#' obj <- mph(structure = c("general", "general"))
#' cdf(obj, matrix(c(0.5, 1), ncol = 2))
setMethod("cdf", c(x = "mph"), function(x,
                                        y,
                                        lower.tail = TRUE) {
  alpha <- x@pars$alpha
  S <- x@pars$S
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
        aux[, i] <- phcdf(y[, i], in_vect, S[[i]], lower.tail)
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
        for (m in 1:n) {
          aux[m, i, j] <- phcdf(y[m, i], in_vect, S[[i]], lower_tail = F)
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

#' Laplace method for multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{mph}.
#' @param r A matrix of real values.
#'
#' @return A vector containing the corresponding Laplace transform evaluations.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- mph(structure = c("general", "general"))
#' laplace(obj, matrix(c(0.5, 1), ncol = 2))
setMethod("laplace", c(x = "mph"), function(x, r) {
  if (methods::is(x, "miph")) {
    warning("Laplace transform of undelying mph structure is provided for miph objects")
  }
  alpha <- x@pars$alpha
  S <- x@pars$S
  d <- length(x@pars$S)

  if (is.matrix(r)) {
    n <- nrow(r)
  }
  if (is.vector(r)) {
    n <- 1
    r <- t(r)
  }

  for (i in 1:d) {
    lim <- max(Re(eigen(S[[i]])$values))
    if (any(r[, i] <= lim)) {
      stop("r should be above the largest real eigenvalue of S")
    }
  }

  res <- numeric(n)

  p <- length(x@pars$alpha)

  for (j in 1:p) {
    in_vect <- rep(0, p)
    in_vect[j] <- 1
    aux <- matrix(NA, n, d)
    for (i in 1:d) {
      for (m in 1:n) {
        aux[m, i] <- ph_laplace(r[m, i], in_vect, S[[i]])
      }
    }
    res <- res + alpha[j] * apply(aux, 1, prod)
  }

  res
})

#' Mgf method for multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{mph}.
#' @param r A matrix of real values.
#'
#' @return A vector containing the corresponding mgf evaluations.
#' @export
#'
#' @examples
#' set.seed(124)
#' obj <- mph(structure = c("general", "general"))
#' mgf(obj, matrix(c(0.5, 0.3), ncol = 2))
setMethod("mgf", c(x = "mph"), function(x, r) {
  if (methods::is(x, "miph")) {
    warning("mgf of undelying mph structure is provided for miph objects")
  }
  S <- x@pars$S
  d <- length(x@pars$S)

  if (is.vector(r)) {
    r <- t(r)
  }

  for (i in 1:d) {
    lim <- -max(Re(eigen(S[[i]])$values))
    if (any(r[, i] > lim)) {
      stop("r should be below the negative largest real eigenvalue of S")
    }
  }

  laplace(x, -r)
})

#' Moment method for multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{mph}.
#' @param k A vector of non-negative integer values.
#'
#' @return The corresponding joint moment evaluation.
#' @export
#'
#' @examples
#' obj <- mph(structure = c("general", "general"))
#' moment(obj, c(2, 1))
setMethod("moment", c(x = "mph"), function(x, k) {
  if (all(k == 0) | any(k < 0)) {
    stop("k should be non-negative and not zero")
  }
  if (any((k %% 1) != 0)) {
    stop("k should be an integer")
  }
  if (methods::is(x, "miph")) {
    warning("moment of undelying mph structure is provided for miph objects")
  }
  alpha <- x@pars$alpha
  S <- x@pars$S
  d <- length(x@pars$S)

  p <- length(x@pars$alpha)
  res <- 0

  for (j in 1:p) {
    in_vect <- rep(0, p)
    in_vect[j] <- 1
    aux <- rep(0, d)
    for (i in 1:d) {
      aux[i] <- factorial(k[i]) * sum(in_vect %*% matrix_power(k[i], solve(-S[[i]])))
    }
    res <- res + alpha[j] * prod(aux)
  }

  res
})

#' Mean method for multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{mph}.
#'
#' @return The mean of the multivariate phase-type distribution.
#' @export
#'
#' @examples
#' obj <- mph(structure = c("general", "general"))
#' mean(obj)
setMethod("mean", c(x = "mph"), function(x) {
  if (methods::is(x, "miph")) {
    warning("mean of undelying mph structure is provided for miph objects")
  }
  d <- length(x@pars$S)
  res <- rep(0, d)
  for (i in 1:d) {
    mom_vect <- rep(0, d)
    mom_vect[i] <- 1
    res[i] <- suppressWarnings(moment(x, mom_vect))
  }

  res
})

#' Var method for multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{mph}.
#'
#' @return The covariance matrix of the multivariate phase-type distribution.
#' @export
#'
#' @examples
#' obj <- mph(structure = c("general", "general"))
#' var(obj)
setMethod("var", c(x = "mph"), function(x) {
  if (methods::is(x, "miph")) {
    warning("covariance matrix of undelying mph structure is provided for miph objects")
  }
  d <- length(x@pars$S)
  res <- matrix(0, d, d)
  for (i in 1:d) {
    mom_vect1 <- rep(0, d)
    mom_vect1[i] <- 1
    for (j in i:d) {
      mom_vect2 <- rep(0, d)
      mom_vect2[j] <- 1
      res[i, j] <- suppressWarnings(moment(x, mom_vect1 + mom_vect2) - moment(x, mom_vect1) * moment(x, mom_vect2))
    }
  }
  res[lower.tri(res)] <- t(res)[lower.tri(res)]
  res
})

#' Cor method for multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{mph}.
#'
#' @return The correlation matrix of the multivariate phase-type distribution.
#' @export
#'
#' @examples
#' obj <- mph(structure = c("general", "general"))
#' cor(obj)
setMethod("cor", c(x = "mph"), function(x) {
  if (methods::is(x, "miph")) {
    warning("correlation matrix of undelying mph structure is provided for miph objects")
  }
  suppressWarnings(stats::cov2cor(var(x)))
})

#' Fit method for mph Class
#'
#' @param x An object of class \linkS4class{mph}.
#' @param y Matrix of data.
#' @param delta Matrix with right-censoring indicators (1 uncensored, 0 right censored).
#' @param stepsEM Number of EM steps to be performed.
#' @param equal_marginals Logical. If `TRUE`, all marginals are fitted to be equal.
#' @param r Sub-sampling parameter, defaults to 1.
#' @param maxit Maximum number of iterations when optimizing g function.
#' @param reltol Relative tolerance when optimizing g function.
#'
#' @export
#'
#' @examples
#' obj <- mph(structure = c("general", "coxian"))
#' data <- sim(obj, 100)
#' fit(x = obj, y = data, stepsEM = 20)
setMethod(
  "fit", c(x = "mph", y = "ANY"),
  function(x, y,
           delta = numeric(0),
           stepsEM = 1000,
           equal_marginals = FALSE,
           r = 1,
           maxit = 100,
           reltol = 1e-8) {
    if (any(y < 0)) {
      stop("data should be positive")
    }
    if (stepsEM <= 0) {
      stop("the number of steps should be positive")
    }
    if (r <= 0 && r > 1) {
      stop("sub-sampling proportion is invalid, please input a r in (0,1]")
    }
    d <- length(x@pars$S)

    is_miph <- methods::is(x, "miph")
    if (is_miph) {
      par_name <- x@gfun$name
      par_g <- x@gfun$pars
      inv_g <- x@gfun$inverse

      opt_fun <- miph_LL
    }

    if (length(delta) == 0) {
      delta <- matrix(1, nrow(y), ncol(y))
    }
    alpha_fit <- x@pars$alpha
    S_fit <- x@pars$S

    if (r < 1) {
      y_full <- y
      delta_full <- delta
    }

    if (!equal_marginals) {
      S_fit <- x@pars$S
      fnn <- EM_step_mph
    } else {
      for (i in 1:ncol(y)) {
        S_fit <- x@pars$S[[1]]
      }
      fnn <- EM_step_mph_0
    }

    options(digits.secs = 4)
    cat(format(Sys.time(), format = "%H:%M:%OS"), ": EM started", sep = "")
    cat("\n", sep = "")
    # EM step
    if (!is_miph) {
      for (k in 1:stepsEM) {
        if (r < 1) {
          index <- sample(1:nrow(y_full), size = floor(r * nrow(y_full)))

          y <- as.matrix(y_full[index, ])
          delta <- as.matrix(delta_full[index, ])
        }

        # EM_step_mPH_rc(alpha_fit, S_fit, y ,delta, h) # C++
        aux <- fnn(alpha_fit, S_fit, y, delta)

        alpha_fit <- aux$alpha
        S_fit <- aux$S
        cat("\r", "iteration:", k,
          ", logLik:", aux$logLik,
          sep = " "
        )
      }
      x@pars$alpha <- alpha_fit # C++
      x@pars$S <- S_fit # C++
      x@fit <- list(
        logLik = sum(log(dens(x, y, delta))),
        nobs = nrow(y)
      )

      if (equal_marginals) {
        ls <- list()
        for (i in 1:ncol(y)) {
          ls[[i]] <- S_fit
        }
        x@pars$S <- ls
      }
    }
    if (is_miph) {
      for (k in 1:stepsEM) {
        # sub-sampling
        if (r < 1) {
          index <- sample(1:nrow(y_full), size = floor(r * nrow(y_full)))

          y <- as.matrix(y_full[index, ])
          delta <- as.matrix(delta_full[index, ])
        }

        # transform to time-homogeneous
        trans <- clone_matrix(y)
        for (i in 1:d) {
          if (x@gfun$name[i] != "gev") {
            trans[, i] <- inv_g[[i]](par_g[[i]], y[, i])
          } else {
            t <- inv_g[[i]](par_g[[i]], y[, i], rep(1, nrow(y)))
            trans[, i] <- t$obs
          }
        }

        aux <- fnn(alpha_fit, S_fit, trans, delta)


        alpha_fit <- aux$alpha
        S_fit <- aux$S

        x@pars$alpha <- alpha_fit
        x@pars$S <- S_fit

        opt <- suppressWarnings(
          stats::optim(
            par = par_g,
            fn = opt_fun,
            x = x,
            obs = y,
            delta = delta,
            hessian = F,
            control = list(
              maxit = maxit,
              reltol = reltol,
              fnscale = -1
            )
          )
        )

        par_g <- as.list(opt$par)

        cat("\r", ", iteration:", k,
          ", logLik:", opt$value,
          sep = " "
        )

        alpha_fit <- aux$alpha
        S_fit <- aux$S
      }
      x@pars$alpha <- alpha_fit
      x@pars$S <- S_fit
      x@gfun$pars <- par_g
      # x <- miph(mph=x, gfun = par_name, gfun_pars= par_g)

      x@fit <- list(
        logLik = sum(log(dens(x, y, delta))),
        nobs = nrow(y)
      )
    }
    cat("\n", sep = "")

    cat("\n", format(Sys.time(), format = "%H:%M:%OS"), ": EM finalized", sep = "")
    cat("\n", sep = "")

    x
  }
)

# multivariate loglikelihood to be optimized
miph_LL <- function(x,
                    obs,
                    delta,
                    gfun_pars) {
  x@gfun$pars <- gfun_pars
  res <- dens(x = x, y = obs, delta = delta)

  sum(log(res))
}

# EM step for mPH class
EM_step_mph <- function(alpha, S_list, y, delta) {
  p <- length(alpha)
  n <- nrow(y)
  d <- ncol(y)
  matrix_integrals <- list()
  for (i in 1:d) {
    s <- -rowSums(S_list[[i]])
    marg <- list()
    for (j in 1:p) {
      e <- rep(0, p)
      e[j] <- 1
      big_mat <- rbind(
        cbind(S_list[[i]], outer(s, e)),
        cbind(matrix(0, p, p), S_list[[i]])
      )

      big_mat_rc <- rbind(
        cbind(S_list[[i]], outer(rep(1, p), e)),
        cbind(matrix(0, p, p), S_list[[i]])
      )

      marg[[j]] <- vapply(
        X = y[, i], FUN = function(yy) {
          matrix_exponential(yy * big_mat)
        },
        FUN.VALUE = matrix(1, 2 * p, 2 * p)
      )
      if (any(delta == 0)) {
        rc <- which(delta[, i] == 0)
        for (m in rc) {
          marg[[j]][, , m] <- matrix_exponential(y[m, i] * big_mat_rc)
        }
      }
    }

    matrix_integrals[[i]] <- marg
  }
  ###
  a_kij <- array(NA, c(n, p, d, p))
  for (k in 1:p) {
    for (j in 1:p) {
      for (i in 1:d) {
        for (m in 1:n) {
          a_kij[m, k, i, j] <- matrix_integrals[[i]][[1]][k, j, m]
        }
      }
    }
  }
  ###
  a_ki <- array(NA, c(n, p, d))
  for (k in 1:p) {
    for (i in 1:d) {
      ti <- -rowSums(S_list[[i]])
      for (m in 1:n) {
        if (delta[m, i] == 1) {
          a_ki[m, k, i] <- sum(ti * a_kij[m, k, i, ])
        } else {
          a_ki[m, k, i] <- sum(a_kij[m, k, i, ])
        }
      }
    }
  }
  ###
  a_k_minus_i <- array(NA, c(n, p, d))
  for (k in 1:p) {
    for (i in 1:d) {
      for (m in 1:n) {
        a_k_minus_i[m, k, i] <- prod(a_ki[m, k, -i])
      }
    }
  }
  ###
  a_k <- array(NA, c(n, p))
  for (k in 1:p) {
    for (m in 1:n) {
      a_k[m, k] <- prod(a_ki[m, k, ])
    }
  }
  ###
  a_tilde_ki <- array(NA, c(n, p, d))
  for (k in 1:p) {
    for (i in 1:d) {
      for (m in 1:n) {
        if (delta[m, i] == 1) {
          a_tilde_ki[m, k, i] <- sum(alpha * a_kij[m, , i, k] * a_k_minus_i[m, , i])
        } else {
          a_tilde_ki[m, k, i] <- 0
        }
      }
    }
  }
  ###
  a <- array(NA, c(n))
  for (m in 1:n) {
    a[m] <- sum(alpha * a_k[m, ])
  }
  ###
  b_skij <- array(NA, c(n, p, p, d, p))
  for (s in 1:p) {
    for (k in 1:p) {
      for (j in 1:p) {
        for (i in 1:d) {
          for (m in 1:n) {
            b_skij[m, s, k, i, j] <- matrix_integrals[[i]][[j]][s, k + p, m]
          }
        }
      }
    }
  }
  ###
  b_ski <- array(NA, c(n, p, p, d))
  for (s in 1:p) {
    for (k in 1:p) {
      for (i in 1:d) {
        for (m in 1:n) {
          b_ski[m, s, k, i] <- sum(alpha * a_k_minus_i[m, , i] * b_skij[m, s, k, i, ])
        }
      }
    }
  }
  ###
  # E STEP
  ###
  EB_k <- numeric(p)
  for (k in 1:p) {
    EB_k[k] <- alpha[k] * sum(a_k[, k] / a)
  }
  EZ_ki <- array(NA, c(p, d))
  for (k in 1:p) {
    for (i in 1:d) {
      EZ_ki[k, i] <- sum(b_ski[, k, k, i] / a)
    }
  }
  EN_ksi <- array(NA, c(p, p, d))
  for (i in 1:d) {
    for (s in 1:p) {
      for (k in 1:p) {
        t_ksi <- S_list[[i]][k, s]
        EN_ksi[k, s, i] <- t_ksi * sum(b_ski[, s, k, i] / a)
      }
    }
  }
  EN_ki <- array(NA, c(p, d))
  for (i in 1:d) {
    for (k in 1:p) {
      t_ki <- -rowSums(S_list[[i]])[k]
      EN_ki[k, i] <- t_ki * sum(a_tilde_ki[, k, i] / a)
    }
  }
  ###
  # M STEP
  ###
  alpha <- numeric(p)
  for (k in 1:p) {
    alpha[k] <- EB_k[k] / n
  }
  S <- array(NA, c(p, p, d))
  for (i in 1:d) {
    for (s in 1:p) {
      for (k in 1:p) {
        S[k, s, i] <- EN_ksi[k, s, i] / EZ_ki[k, i]
      }
    }
  }
  s <- array(NA, c(p, d))
  for (i in 1:d) {
    for (k in 1:p) {
      s[k, i] <- EN_ki[k, i] / EZ_ki[k, i]
    }
  }
  for (k in 1:p) {
    for (i in 1:d) {
      S[k, k, i] <- -sum(S[k, -k, i]) - s[k, i]
    }
  }
  ll <- list()
  for (i in 1:d) {
    ll[[i]] <- S[, , i]
  }
  list(alpha = alpha, S = ll, logLik = sum(log(a)))
}

EM_step_mph_0 <- function(alpha, S, y, delta) {
  p <- length(alpha)
  n <- nrow(y)
  d <- ncol(y)
  matrix_integrals <- list()
  s <- -rowSums(S)
  S_list <- NULL # Added to avoid a Note in the check - Not sure
  for (i in 1:d) {
    marg <- list()
    for (j in 1:p) {
      e <- rep(0, p)
      e[j] <- 1
      big_mat <- rbind(
        cbind(S, outer(s, e)),
        cbind(matrix(0, p, p), S)
      )
      big_mat_rc <- rbind(
        cbind(S_list[[i]], outer(rep(1, p), e)),
        cbind(matrix(0, p, p), S_list[[i]])
      )

      marg[[j]] <- vapply(
        X = y[, i], FUN = function(yy) {
          matrix_exponential(yy * big_mat)
        },
        FUN.VALUE = matrix(1, 2 * p, 2 * p)
      )
      if (any(delta == 0)) {
        rc <- which(delta[, i] == 0)
        for (m in rc) {
          marg[[j]][, , m] <- matrix_exponential(y[m, i] * big_mat_rc)
        }
      }
    }
    matrix_integrals[[i]] <- marg
  }
  ###
  a_kij <- array(NA, c(n, p, d, p))
  for (k in 1:p) {
    for (j in 1:p) {
      for (i in 1:d) {
        for (m in 1:n) {
          a_kij[m, k, i, j] <- matrix_integrals[[i]][[1]][k, j, m]
        }
      }
    }
  }
  ###
  a_ki <- array(NA, c(n, p, d))
  for (k in 1:p) {
    for (i in 1:d) {
      ti <- -rowSums(S_list[[i]])
      for (m in 1:n) {
        if (delta[m, i] == 1) {
          a_ki[m, k, i] <- sum(ti * a_kij[m, k, i, ])
        } else {
          a_ki[m, k, i] <- sum(a_kij[m, k, i, ])
        }
      }
    }
  }
  ###
  a_k_minus_i <- array(NA, c(n, p, d))
  for (k in 1:p) {
    for (i in 1:d) {
      for (m in 1:n) {
        a_k_minus_i[m, k, i] <- prod(a_ki[m, k, -i])
      }
    }
  }
  ###
  a_k <- array(NA, c(n, p))
  for (k in 1:p) {
    for (m in 1:n) {
      a_k[m, k] <- prod(a_ki[m, k, ])
    }
  }
  ###
  a_tilde_ki <- array(NA, c(n, p, d))
  for (k in 1:p) {
    for (i in 1:d) {
      for (m in 1:n) {
        if (delta[m, i] == 1) {
          a_tilde_ki[m, k, i] <- sum(alpha * a_kij[m, , i, k] * a_k_minus_i[m, , i])
        } else {
          a_tilde_ki[m, k, i] <- 0
        }
      }
    }
  }
  ###
  a <- array(NA, c(n))
  for (m in 1:n) {
    a[m] <- sum(alpha * a_k[m, ])
  }
  ###
  b_skij <- array(NA, c(n, p, p, d, p))
  for (s in 1:p) {
    for (k in 1:p) {
      for (j in 1:p) {
        for (i in 1:d) {
          for (m in 1:n) {
            b_skij[m, s, k, i, j] <- matrix_integrals[[i]][[j]][s, k + p, m]
          }
        }
      }
    }
  }
  ###
  b_ski <- array(NA, c(n, p, p, d))
  for (s in 1:p) {
    for (k in 1:p) {
      for (i in 1:d) {
        for (m in 1:n) {
          b_ski[m, s, k, i] <- sum(alpha * a_k_minus_i[m, , i] * b_skij[m, s, k, i, ])
        }
      }
    }
  }
  ###
  # E STEP
  ###
  EB_k <- numeric(p)
  for (k in 1:p) {
    EB_k[k] <- alpha[k] * sum(a_k[, k] / a)
  }
  EZ_ki <- array(NA, c(p, d))
  for (k in 1:p) {
    for (i in 1:d) {
      EZ_ki[k, i] <- sum(b_ski[, k, k, i] / a)
    }
  }
  EZ_k <- rowSums(EZ_ki, dims = 1)
  EN_ksi <- array(NA, c(p, p, d))
  for (i in 1:d) {
    for (s in 1:p) {
      for (k in 1:p) {
        t_ksi <- S[k, s]
        EN_ksi[k, s, i] <- t_ksi * sum(b_ski[, s, k, i] / a)
      }
    }
  }
  EN_ks <- rowSums(EN_ksi, dims = 2)
  EN_ki <- array(NA, c(p, d))
  for (i in 1:d) {
    for (k in 1:p) {
      t_ki <- -rowSums(S)[k]
      EN_ki[k, i] <- t_ki * sum(a_tilde_ki[, k, i] / a)
    }
  }
  EN_k <- rowSums(EN_ki, dims = 1)
  ###
  # M STEP
  ###
  alpha <- numeric(p)
  for (k in 1:p) {
    alpha[k] <- EB_k[k] / n
  }
  S <- matrix(NA, p, p)
  for (s in 1:p) {
    for (k in 1:p) {
      S[k, s] <- EN_ks[k, s] / EZ_k[k]
    }
  }
  s <- numeric(p)
  for (k in 1:p) {
    s[k] <- EN_k[k] / EZ_k[k]
  }
  for (k in 1:p) {
    S[k, k] <- -sum(S[k, -k]) - s[k]
  }
  list(alpha = alpha, S = S, logLik = sum(log(a)))
}
