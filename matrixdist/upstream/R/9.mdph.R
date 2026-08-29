#' Multivariate discrete phase-type distributions
#'
#' Class of objects for multivariate discrete phase-type distributions.
#'
#' @slot name Name of the discrete phase-type distribution.
#' @slot pars A list comprising of the parameters.
#' @slot fit A list containing estimation information.
#'
#' @return Class object.
#' @export
#'
setClass("mdph",
  slots = list(
    name = "character",
    pars = "list",
    fit = "list"
  )
)

#' Constructor function for multivariate discrete phase-type distributions
#'
#' @param alpha A probability vector.
#' @param S A list of sub-transition matrices.
#' @param structure A vector of valid ph structures.
#' @param dimension The dimension of the dph structure (if provided).
#' @param variables The dimension of the multivariate discrete phase-type.
#'
#' @return An object of class \linkS4class{mdph}.
#' @export
#'
#' @examples
#' mdph(structure = c("general", "general"), dimension = 5)
mdph <- function(alpha = NULL, S = NULL, structure = NULL, dimension = 3, variables = NULL) {
  if (any(is.null(alpha)) & any(is.null(S)) & is.null(structure)) {
    stop("input a vector and matrix, or a structure")
  }
  if (is.null(variables)) {
    variables <- length(structure)
  }
  if (!any(is.null(structure))) {
    rs <- random_structure(dimension, structure = structure[1])
    alpha <- rs[[1]]
    Sa <- rs[[2]]
    a <- max_diagonal(Sa * (-1)) * (1 + stats::runif(1))
    S <- list()
    S[[1]] <- (a * diag(dimension) + Sa) / a
    for (i in 2:variables) {
      Sa <- random_structure(dimension, structure = structure[i])[[2]]
      a <- max_diagonal(Sa * (-1)) * (1 + stats::runif(1))
      S[[i]] <- (a * diag(dimension) + Sa) / a
    }
    name <- structure
  } else {
    name <- "custom"
  }
  methods::new("mdph",
    name = paste(name, " mdph(", length(alpha), ")", sep = " "),
    pars = list(alpha = alpha, S = S)
  )
}

#' Show method for multivariate discrete phase-type distributions
#'
#' @param object An object of class \linkS4class{mdph}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "mdph", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@pars)
  cat("number of variables: ", length(object@pars$S), "\n", sep = "")
})

#' Coef method for mdph class
#'
#' @param object An object of class \linkS4class{mdph}.
#'
#' @return Parameters of multivariate discrete phase-type model.
#' @export
#'
#' @examples
#' obj <- mdph(structure = c("general", "general"))
#' coef(obj)
setMethod("coef", c(object = "mdph"), function(object) {
  object@pars
})

#' Simulation method for multivariate discrete phase-type distributions
#'
#' @param x An object of class \linkS4class{mdph}.
#' @param n Length of realization.
#' @param equal_marginals Non-negative integer. If positive, it specifies
#' the number of marginals to simulate from, all from the first matrix.
#'
#' @return A realization of a multivariate discrete phase-type distribution.
#' @export
#'
#' @examples
#' obj <- mdph(structure = c("general", "general"))
#' sim(obj, 100)
setMethod("sim", c(x = "mdph"), function(x, n = 1000, equal_marginals = 0) {
  p <- length(x@pars$alpha)

  if (equal_marginals == 0) {
    d <- length(x@pars$S)

    trans_mat <- list()
    for (j in 1:d) {
      exit_vec <- 1 - rowSums(x@pars$S[[j]])
      t_mat <- cbind(x@pars$S[[j]], exit_vec)
      aux_vec <- rep(0, p + 1)
      aux_vec[p + 1] <- 1
      trans_mat[[j]] <- rbind(t_mat, aux_vec)
    }

    states <- 1:p
    result <- matrix(NA, n, d)
    for (i in 1:n) {
      state <- sample(states, 1, prob = x@pars$alpha)
      in_vect <- rep(0, p)
      in_vect[state] <- 1
      for (j in 1:d) {
        result[i, j] <- rdphasetype(1, in_vect, trans_mat[[j]])
      }
    }
  } else {
    d <- equal_marginals

    exit_vec <- 1 - rowSums(x@pars$S[[1]])
    trans_mat <- cbind(x@pars$S[[1]], exit_vec)
    aux_vec <- rep(0, p + 1)
    aux_vec[p + 1] <- 1
    trans_mat <- rbind(trans_mat, aux_vec)

    states <- 1:p
    result <- matrix(NA, n, d)
    for (i in 1:n) {
      state <- sample(states, 1, prob = x@pars$alpha)
      in_vect <- rep(0, p)
      in_vect[state] <- 1
      for (j in 1:d) {
        result[i, j] <- rdphasetype(1, in_vect, trans_mat)
      }
    }
  }
  result
})

#' Marginal method for mdph class
#'
#' @param x An object of class \linkS4class{mdph}.
#' @param mar Indicator of which marginal.
#' @return An object of the of class \linkS4class{dph}.
#' @export
#'
#' @examples
#' obj <- mdph(structure = c("general", "general"))
#' marginal(obj, 1)
setMethod("marginal", c(x = "mdph"), function(x, mar = 1) {
  if (!(mar %in% 1:length(x@pars$S))) {
    stop("maringal provided not available")
  }
  dph(alpha = x@pars$alpha, S = x@pars$S[[mar]])
})

#' Moment method for multivariate discrete phase-type distributions
#'
#' @param x An object of class \linkS4class{mdph}.
#' @param k A vector of positive integer values.
#'
#' @return The corresponding joint factorial moment evaluation.
#' @export
#'
#' @examples
#' obj <- mdph(structure = c("general", "general"))
#' moment(obj, c(2, 1))
setMethod("moment", c(x = "mdph"), function(x, k) {
  if (any(k <= 0)) {
    stop("k should have positive entries")
  }
  if (any((k %% 1) != 0)) {
    stop("k should be an integer")
  }

  d <- length(x@pars$S)

  if (length(k) != d) {
    stop("vector k of wrong dimension")
  }

  alpha <- x@pars$alpha
  S <- x@pars$S

  p <- length(x@pars$alpha)
  res <- 0

  for (j in 1:p) {
    in_vect <- rep(0, p)
    in_vect[j] <- 1
    aux <- rep(0, d)
    for (i in 1:d) {
      m1 <- matrix_power(k[i] - 1, S[[i]])
      m2 <- matrix_power(k[i], solve(diag(p) - S[[i]]))
      aux[i] <- factorial(k[i]) * sum(in_vect %*% m1 %*% m2)
    }
    res <- res + alpha[j] * prod(aux)
  }

  res
})

#' Mean method for multivariate discrete phase-type distributions
#'
#' @param x An object of class \linkS4class{mdph}.
#'
#' @return The mean of the multivariate discrete phase-type distribution.
#' @export
#'
#' @examples
#' obj <- mdph(structure = c("general", "general"))
#' mean(obj)
setMethod("mean", c(x = "mdph"), function(x) {
  d <- length(x@pars$S)
  res <- rep(0, d)
  for (i in 1:d) {
    res[i] <- mean(marginal(x, i))
  }
  res
})

#' Var method for multivariate discrete phase-type distributions
#'
#' @param x An object of class \linkS4class{mdph}.
#'
#' @return The covariance matrix of the multivariate discrete phase-type distribution.
#' @export
#'
#' @examples
#' obj <- mdph(structure = c("general", "general"))
#' var(obj)
setMethod("var", c(x = "mdph"), function(x) {
  d <- length(x@pars$S)
  alpha <- x@pars$alpha
  S <- x@pars$S
  p <- length(x@pars$alpha)

  res <- matrix(0, d, d)
  for (i in 1:d) {
    mar1 <- marginal(x, i)
    for (j in i:d) {
      if (j == i) {
        res[i, j] <- var(mar1)
      } else {
        mar2 <- marginal(x, j)
        cross <- 0
        for (s in 1:p) {
          in_vect <- rep(0, p)
          in_vect[s] <- 1
          aux <- rep(0, 2)

          m1 <- solve(diag(p) - S[[i]])
          aux[1] <- sum(in_vect %*% m1)

          m1 <- solve(diag(p) - S[[j]])
          aux[2] <- sum(in_vect %*% m1)

          cross <- cross + alpha[s] * prod(aux)
        }
        res[i, j] <- cross - mean(mar1) * mean(mar2)
      }
    }
  }
  res[lower.tri(res)] <- t(res)[lower.tri(res)]
  res
})

#' Cor method for multivariate discrete phase-type distributions
#'
#' @param x An object of class \linkS4class{mdph}.
#'
#' @return The correlation matrix of the multivariate discrete phase-type distribution.
#' @export
#'
#' @examples
#' obj <- mdph(structure = c("general", "general"))
#' cor(obj)
setMethod("cor", c(x = "mdph"), function(x) {
  stats::cov2cor(var(x))
})

#' Pgf method for multivariate discrete phase-type distributions
#'
#' @param x An object of class \linkS4class{mdph}.
#' @param z A matrix of real values.
#'
#' @return A vector containing the corresponding pgf evaluations.
#' @export
#'
#' @examples
#' obj <- mdph(structure = c("general", "general"))
#' pgf(obj, matrix(c(0.5, 1), ncol = 2))
setMethod("pgf", c(x = "mdph"), function(x, z) {
  alpha <- x@pars$alpha
  S <- x@pars$S
  d <- length(x@pars$S)

  if (is.matrix(z)) {
    n <- nrow(z)
  }
  if (is.vector(z)) {
    n <- 1
    z <- t(z)
  }

  if (any(abs(z) > 1)) {
    stop("z should between -1 and 1")
  }

  res <- numeric(n)

  p <- length(x@pars$alpha)

  for (j in 1:p) {
    in_vect <- rep(0, p)
    in_vect[j] <- 1
    aux <- matrix(NA, n, d)
    for (i in 1:d) {
      for (m in 1:n) {
        aux[m, i] <- dph_pgf(z[m, i], in_vect, S[[i]])
      }
    }
    res <- res + alpha[j] * apply(aux, 1, prod)
  }

  res
})

#' Density method for multivariate discrete phase-type distributions
#'
#' @param x An object of class \linkS4class{mdph}.
#' @param y A matrix of locations.
#'
#' @return A vector containing the joint density evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- mdph(structure = c("general", "general"))
#' dens(obj, matrix(c(1, 1), ncol = 2))
setMethod("dens", c(x = "mdph"), function(x, y) {
  if (is.vector(y)) {
    y <- t(y)
  }
  mdphdensity(y, x@pars$alpha, x@pars$S)
})

#' Fit method for mdph Class
#'
#' @param x An object of class \linkS4class{mdph}.
#' @param y A matrix with the data.
#' @param weight Vector of weights.
#' @param stepsEM Number of EM steps to be performed.
#' @param every Number of iterations between likelihood display updates.
#'
#' @return An object of class \linkS4class{mdph}.
#'
#' @export
#'
#' @examples
#' obj <- mdph(structure = c("general", "general"))
#' data <- sim(obj, n = 100)
#' fit(obj, data, stepsEM = 100, every = 50)
setMethod(
  "fit", c(x = "mdph"),
  function(x,
           y,
           weight = numeric(0),
           stepsEM = 1000,
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

    mdph_par <- x@pars
    alpha_fit <- clone_vector(mdph_par$alpha)
    S_fit <- list()
    for (j in 1:length(y[1, ])) {
      S_fit[[j]] <- clone_matrix(mdph_par$S[[j]])
    }

    options(digits.secs = 4)
    cat(format(Sys.time(), format = "%H:%M:%OS"), ": EM started", sep = "")
    cat("\n", sep = "")

    for (k in 1:stepsEM) {
      EMstep_mdph(alpha_fit, S_fit, y, weight)
      if (k %% every == 0) {
        cat("\r", "iteration:", k,
          ", logLik:", logLikelihoodmDPH(alpha_fit, S_fit, y, weight),
          sep = " "
        )
      }
    }

    x@pars$alpha <- alpha_fit
    x@pars$S <- S_fit

    cat("\n", format(Sys.time(), format = "%H:%M:%OS"), ": EM finalized", sep = "")
    cat("\n", sep = "")

    x
  }
)

#' MoE method for mdph Class
#'
#' @param x An object of class \linkS4class{mdph}.
#' @param formula A regression formula.
#' @param y A matrix of observations.
#' @param data A data frame of covariates.
#' @param alpha_vecs Matrix of initial probabilities.
#' @param weight Vector of weights.
#' @param stepsEM Number of EM steps to be performed.
#' @param every Number of iterations between likelihood display updates.
#' @param rand_init Random initiation in the R-step.
#' @param maxWts Maximal number of weights in the nnet function.
#'
#' @return An object of class \linkS4class{sph}.
#'
#' @export
#'
#' @examples
#' x <- mdph(structure = c("general", "general"))
#' n <- 100
#' responses <- cbind(rpois(n, 3) + 1, rbinom(n, 5, 0.5))
#' covariates <- data.frame(age = sample(18:65, n, replace = TRUE) / 100, income = runif(n, 0, 0.99))
#' f <- responses ~ age + income
#' MoE(x = x, formula = f, y = responses, data = covariates, stepsEM = 20)
setMethod(
  "MoE", c(x = "mdph"),
  function(x,
           formula,
           y,
           data,
           alpha_vecs = NULL,
           weight = numeric(0),
           stepsEM = 1000,
           every = 10,
           rand_init = TRUE,
           maxWts = 1000) {
    p <- length(x@pars$alpha)
    frame <- stats::model.frame(formula, data = data)
    n <- nrow(frame)
    d <- ncol(frame) - 1
    if (is.null(alpha_vecs)) alpha_vecs <- matrix(x@pars$alpha, ncol = p, nrow = n, byrow = TRUE)
    if (length(weight) == 0) weight <- rep(1, n)

    S_fit <- list()
    for (j in 1:length(y[1, ])) {
      S_fit[[j]] <- clone_matrix(x@pars$S[[j]])
    }

    c <- c()
    for (i in 1:p) c <- c(c, rep(i, n))
    extended_x <- matrix(t(as.matrix(frame[, -1])), nrow = n * p, ncol = d, byrow = TRUE)
    dm <- data.frame(Class = c, extended_x)
    names(dm)[-1] <- names(frame)[-1]
    ndm <- data.frame(dm[dm$Class == 1, -1])
    names(ndm) <- names(dm)[-1]
    for (k in 1:stepsEM) {
      B_matrix <- EMstep_mdph_MoE(alpha_vecs, S_fit, y, weight)
      wt <- reshape2::melt(B_matrix)[, 3]
      wt[wt < 1e-22] <- wt[wt < 1e-22] + 1e-22
      if (k == 1 | rand_init == TRUE) {
        multinom_model <- nnet::multinom(Class ~ ., data = dm, weights = wt, trace = F, MaxNWts = maxWts)
      } else {
        multinom_model <- nnet::multinom(Class ~ ., data = dm, weights = wt, trace = F, Wts = multinom_model$wts, MaxNWts = maxWts)
      }
      alpha_vecs <- stats::predict(multinom_model, type = "probs", newdata = ndm)
      if (k %% every == 0) {
        ll <- logLikelihoodmDPH_MoE(alpha_vecs, S_fit, y, weight)
        cat("\r", "iteration:", k, ", logLik:", ll, sep = " ")
      }
    }
    cat("\n", sep = "")
    list(alpha = alpha_vecs, S = S_fit, mm = multinom_model)
  }
)
