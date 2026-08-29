#' Multivariate phase-type distributions obtained by transformation via rewards
#'
#' Class of objects for multivariate phase type distributions.
#'
#' @slot name Name of the phase type distribution.
#' @slot pars A list comprising of the parameters.
#'
#' @export
#'
setClass("MPHstar",
  slots = list(
    name = "character",
    pars = "list"
  )
)

#' Constructor function for multivariate phase-type distributions (MPH* class)
#'
#' @param alpha A probability vector.
#' @param S A sub-intensity matrix.
#' @param R A compatible (non-negative) reward matrix.
#' @param structure A valid ph structure.
#' @param dimension The dimension of the ph structure (if provided).
#' @param variables The number of desired marginals.
#'
#' @return An object of class \linkS4class{MPHstar}.
#' @export
#'
#' @examples
#' MPHstar(structure = "general", dimension = 4, variables = 3)
MPHstar <- function(alpha = NULL,
                    S = NULL,
                    structure = NULL,
                    dimension = 3,
                    R = NULL,
                    variables = 2) {
  if (any(is.null(alpha)) & any(is.null(S)) & is.null(structure)) {
    stop("input a vector and matrix, or a structure")
  }

  if (!all(is.null(R))) {
    if (!all(R %*% rep(1, ncol(R)) == 1)) stop("Rows of the reward matrix must sum to 1")
    Rname <- "custom"
  }

  if (!is.null(structure)) {
    rs <- random_structure(p = dimension, structure = structure, scale_factor = 1)
    alpha <- rs[[1]]
    S <- rs[[2]]
    name <- structure
  } else {
    if (dim(S)[1] != dim(S)[2]) {
      stop("matrix S should be square")
    }
    if (length(alpha) != dim(S)[1]) {
      stop("incompatible dimensions")
    }
    name <- "custom"
  }

  if (all(is.null(R))) {
    R <- random_reward(dim(S)[1], variables)
    Rname <- "random"
  } else {
    if (dim(R)[1] != dim(S)[1]) {
      stop("The reward matrix R has wrong dimensions: number of rows does not match ph dimension")
      Rname <- "custom"
    }
  }

  methods::new("MPHstar",
    name = paste("mph based on a ", name, " ph and a ", Rname, " reward matrix", sep = ""),
    pars = list(alpha = alpha, S = S, R = R)
  )
}

#' Show method for multivariate phase-type distributions
#'
#' @param object An object of class \linkS4class{MPHstar}.
#' @export
#'
setMethod("show", "MPHstar", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@pars)
  cat("number of marginal distributions: ", dim(object@pars$R)[2], "\n", sep = "")
})

#' Simulation method for multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{MPHstar}.
#' @param n Desired sample size for each marginal.
#'
#' @return A matrix of sample data for each marginal.
#' @export
#'
#' @examples
#' obj <- MPHstar(structure = "general")
#' sim(obj, 100)
setMethod("sim", c(x = "MPHstar"), function(x, n = 1000) {
  rMPHstar(n, x@pars$alpha, x@pars$S, x@pars$R)
})

#' Marginal method for MPHstar class
#'
#' @param x An object of class \linkS4class{MPHstar}.
#' @param mar Indicator of which marginal.
#' @return An object of the of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' obj <- MPHstar(structure = "general")
#' marginal(obj, 1)
setMethod("marginal", c(x = "MPHstar"), function(x, mar = 1) {
  if (!(mar %in% 1:ncol(x@pars$R))) {
    stop("maringal provided not available")
  }
  mar_par <- tvr_ph(x@pars$alpha, x@pars$S, x@pars$R[, mar])
  ph(alpha = mar_par[[1]], S = mar_par[[2]])
})

#' Linear combination method for MPHstar class
#'
#' @param x An object of class \linkS4class{MPHstar}.
#' @param w A vector with non-negative entries.
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' obj <- MPHstar(structure = "general")
#' linCom(obj, c(1, 0))
setMethod("linCom", c(x = "MPHstar"), function(x, w) {
  if (length(w) != ncol(x@pars$R)) {
    stop("vector of wrong dimension")
  }
  if (any(w < 0)) {
    stop("vector with negative entries")
  }
  if (all(w == 0)) {
    stop("vector with all entries zero")
  }
  L <- linear_combination(w, x@pars$alpha, x@pars$S, x@pars$R)
  ph(alpha = L$alpha, S = L$S)
})

#' Mean method for MPHstar class
#'
#' @param x An object of class \linkS4class{MPHstar}.
#'
#' @return The mean of MPHstar distribution.
#' @export
#'
#' @examples
#' obj <- MPHstar(structure = "general")
#' mean(obj)
setMethod("mean", c(x = "MPHstar"), function(x) {
  d <- ncol(x@pars$R)
  res <- rep(0, d)
  for (i in 1:d) {
    res[i] <- mean(marginal(x, i))
  }
  res
})

#' Var method for MPHstar class
#'
#' @param x An object of class \linkS4class{MPHstar}.
#'
#' @return The covariance matrix of the MPHstar distribution.
#' @export
#'
#' @examples
#' obj <- MPHstar(structure = "general")
#' var(obj)
setMethod("var", c(x = "MPHstar"), function(x) {
  alpha <- x@pars$alpha
  U <- solve(-x@pars$S)
  R <- x@pars$R
  d <- ncol(R)

  res <- matrix(0, d, d)
  for (i in 1:d) {
    mar1 <- marginal(x, i)
    for (j in i:d) {
      if (j == i) {
        res[i, j] <- var(mar1)
      } else {
        mar2 <- marginal(x, j)
        cross <- alpha %*% U %*% diag(R[, i]) %*% U %*% R[, j] + alpha %*% U %*% diag(R[, j]) %*% U %*% R[, i]
        res[i, j] <- cross - mean(mar1) * mean(mar2)
      }
    }
  }
  res[lower.tri(res)] <- t(res)[lower.tri(res)]
  res
})

#' Cor method for MPHstar class
#'
#' @param x An object of class \linkS4class{MPHstar}.
#'
#' @return The correlation matrix of the MPHstar distribution.
#' @export
#'
#' @examples
#' obj <- MPHstar(structure = "general")
#' cor(obj)
setMethod("cor", c(x = "MPHstar"), function(x) {
  stats::cov2cor(var(x))
})

#' Find weight of observations
#'
#' @param x A vector of observations from which we want to know their weights.
#'
#' @return A matrix with unique observations as first column and associated weights for second column.
#' @export
#'
find_weight <- function(x) {
  cum_weight <- numeric(0)
  x <- sort(x)
  unique <- unique(x)
  for (i in unique) {
    if (i > 0) {
      cum_weight <- c(cum_weight, length(which(x == i)))
    }
  }
  cbind(unique, cum_weight, deparse.level = 0)
}

#' Prepare data for the MPHstar_EMstep_UNI
#'
#' @param y A matrix with marginal observations, each column corresponds to a marginal.
#' @param w A matrix of weights, each column corresponds to a marginal.
#'
#' @return For summed and marginal observations we have a list with matrices of
#'  unique observations and their associated weights, separated by uncensored and right-censored data.
#' @export
#'
MPHstar_data_aggregation <- function(y, w = numeric(0)) {
  mat <- list()
  sumData <- rowSums(y)
  mat[[1]] <- find_weight(sumData)

  for (i in 1:ncol(y)) {
    m <- 1 + i
    m_y <- y[, i]
    if (length(w) == 0) {
      mat[[m]] <- find_weight(m_y)
    }
    if (length(w) > 0) {
      m_w <- w[, i]

      mat[[m]] <- cbind(m_y, m_w)
    }
  }
  mat
}

#' Fit method for mph class
#'
#' @param x An object of class \linkS4class{MPHstar}.
#' @param y A matrix of marginal data.
#' @param weight A matrix of marginal weights.
#' @param stepsEM The number of EM steps to be performed, defaults to 1000.
#' @param uni_epsilon The epsilon parameter for the uniformization method, defaults to 1e-4.
#' @param zero_tol The smallest value that a reward can take (to avoid numerical instability), defaults to 1e-4.
#' @param every The number of iterations between likelihood display updates. The originating distribution is used, given that there is no explicit density.
#' @param plot Boolean that determines if the plot of the loglikelihood evolution is plotted, defaults to False.
#' @param r The sub-sampling proportion for stochastic EM, defaults to 1.
#' @param replace Boolean that determines if sub-sampling is done with replacement or not, defaults to False.
#'
#' @return An object of class \linkS4class{MPHstar}.
#'
#' @importFrom grDevices dev.off
#' @importFrom graphics hist legend lines
#' @importFrom utils head
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- MPHstar(structure = "general")
#' data <- sim(obj, 100)
#' fit(obj, data, stepsEM = 20)
setMethod(
  "fit", c(x = "MPHstar", y = "ANY"),
  function(x,
           y,
           weight = numeric(0),
           stepsEM = 1000,
           uni_epsilon = 1e-4,
           zero_tol = 1e-4,
           every = 100,
           plot = F,
           r = 1,
           replace = F) {
    if (!is.matrix(y)) {
      stop("Please use a matrix with marginal observations")
    }
    if (!all(y > 0)) {
      stop("data should be positive")
    }
    if (!all(weight >= 0)) {
      stop("weights should be non-negative")
    }
    if (!all(x@pars$R >= 0)) {
      stop("rewards should be non-negative")
    }
    if (all(x@pars$R == 0)) {
      stop("A reward matrix of zeros results in an atom at zero with probability 100%")
    }
    if (r <= 0 && r > 1) {
      stop("sub-sampling proportion is invalid, please input a r in (0,1]")
    }

    is_mph <- methods::is(x, "MPHstar")

    A <- MPHstar_data_aggregation(y, weight)
    B <- A
    C <- A

    ph_par <- x@pars
    alpha_fit <- clone_vector(ph_par$alpha)
    S_fit <- clone_matrix(ph_par$S)
    R_fit <- clone_matrix(ph_par$R)

    log_check <- rep(0, stepsEM / every)

    if (is_mph) {
      options(digits.secs = 4)
      cat(format(Sys.time(), format = "%H:%M:%OS"), ": EM started", sep = "")
      cat("\n", sep = "")

      for (k in 1:stepsEM) {
        if (r < 1) {
          r_weight <- numeric(0)

          indices <- sample(1:nrow(y), size = floor(r * nrow(y)), replace = replace)
          r_y <- y[indices, ]

          if (length(weight) > 0) {
            r_weight <- weight[indices, ]
          }

          B <- MPHstar_data_aggregation(r_y, r_weight)
        }

        epsilon <- if (!is.na(uni_epsilon)) {
          uni_epsilon
        } else {
          1e-4
        }

        MPHstar_EMstep_UNI(epsilon, zero_tol, alpha_fit, S_fit, R_fit, B) # performs a EM step, changes alpha, S and R

        log_lik <- logLikelihoodPH_UNI(epsilon, alpha_fit, S_fit, C[[1]][, 1], C[[1]][, 2], numeric(0), numeric(0))

        log_check[k] <- log_lik

        if (k %% every == 0) {
          cat("\r", "iteration:", k,
            ", logLik:", log_lik,
            sep = " "
          )
        }
      }
    }
    if (plot == T) {
      plot(log_check, type = "l", main = "Evolution of incomplete loglikelihood", ylab = "loglik", xlab = "step EM")
    }

    cat("\n", sep = "")
    x@pars$alpha <- alpha_fit
    x@pars$S <- S_fit
    x@pars$R <- R_fit

    cat("\n", format(Sys.time(), format = "%H:%M:%OS"), ": EM finalized", sep = "")
    cat("\n", sep = "")

    x
  }
)
