#' Discrete phase-type distributions
#'
#' Class of objects for discrete phase-type distributions.
#'
#' @slot name Name of the discrete phase-type distribution.
#' @slot pars A list comprising of the parameters.
#' @slot fit A list containing estimation information.
#'
#' @return Class object.
#' @export
#'
setClass("dph",
  slots = list(
    name = "character",
    pars = "list",
    fit = "list"
  )
)

#' Constructor function for discrete phase-type distributions
#'
#' @param alpha A probability vector.
#' @param S A sub-transition matrix.
#' @param structure A valid dph structure: `"general"`, `"coxian"`,
#'  `"hyperexponential"`, "gcoxian", or `"gerlang"`.
#' @param dimension The dimension of the dph structure (if structure is provided).
#'
#' @return An object of class \linkS4class{dph}.
#' @export
#'
#' @examples
#' dph(structure = "general", dim = 5)
#' dph(alpha = c(0.5, 0.5), S = matrix(c(0.1, 0.5, 0.5, 0.2), 2, 2))
dph <- function(alpha = NULL, S = NULL, structure = NULL, dimension = 3) {
  if (any(is.null(alpha)) & any(is.null(S)) & is.null(structure)) {
    stop("input a vector and matrix, or a structure")
  }
  if (!is.null(structure)) {
    rs <- random_structure(dimension, structure = structure)
    alpha <- rs[[1]]
    Sa <- rs[[2]]
    a <- max_diagonal(Sa * (-1)) * (1 + stats::runif(1))
    S <- (a * diag(dimension) + Sa) / a
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
  methods::new("dph",
    name = paste(name, " dph(", length(alpha), ")", sep = ""),
    pars = list(alpha = alpha, S = S)
  )
}

#' Sum method for discrete phase-type distributions
#'
#' @param e1 An object of class \linkS4class{dph}.
#' @param e2 An object of class \linkS4class{dph}.
#'
#' @return An object of class \linkS4class{dph}.
#' @export
#'
#' @examples
#' dph1 <- dph(structure = "general", dimension = 3)
#' dph2 <- dph(structure = "general", dimension = 5)
#' dph_sum <- dph1 + dph2
#' dph_sum
setMethod(
  "+", signature(e1 = "dph", e2 = "dph"),
  function(e1, e2) {
    L <- sum_dph(e1@pars$alpha, e1@pars$S, e2@pars$alpha, e2@pars$S)
    dph(alpha = L$alpha, S = L$S)
  }
)

#' Minimum method for discrete phase-type distributions
#'
#' @param x1 An object of class \linkS4class{dph}.
#' @param x2 An object of class \linkS4class{dph}.
#'
#' @return An object of class \linkS4class{dph}.
#' @export
#'
#' @examples
#' dph1 <- dph(structure = "general", dimension = 3)
#' dph2 <- dph(structure = "general", dimension = 5)
#' dph_min <- minimum(dph1, dph2)
#' dph_min
setMethod(
  "minimum", signature(x1 = "dph", x2 = "dph"),
  function(x1, x2) {
    alpha <- kronecker(x1@pars$alpha, x2@pars$alpha)
    S <- kronecker(x1@pars$S, x2@pars$S)
    dph(alpha = alpha, S = S)
  }
)

#' Maximum method for discrete phase-type distributions
#'
#' @param x1 An object of class \linkS4class{dph}.
#' @param x2 An object of class \linkS4class{dph}.
#'
#' @return An object of class \linkS4class{dph}.
#' @export
#'
#' @examples
#' dph1 <- dph(structure = "general", dimension = 3)
#' dph2 <- dph(structure = "general", dimension = 5)
#' dph_max <- maximum(dph1, dph2)
#' dph_max
setMethod(
  "maximum", signature(x1 = "dph", x2 = "dph"),
  function(x1, x2) {
    n1 <- length(x1@pars$alpha)
    n2 <- length(x2@pars$alpha)
    alpha <- c(kronecker(x1@pars$alpha, x2@pars$alpha), rep(0, n1 + n2))
    S1 <- rbind(kronecker(x1@pars$S, x2@pars$S), matrix(0, n1 + n2, n1 * n2))
    S2 <- rbind(kronecker(x1@pars$S, 1 - rowSums(x2@pars$S)), x1@pars$S, matrix(0, n2, n1))
    S3 <- rbind(kronecker(1 - rowSums(x1@pars$S), x2@pars$S), matrix(0, n1, n2), x2@pars$S)
    dph(alpha = alpha, S = cbind(S1, S2, S3))
  }
)

#' Mixture method for phase-type distributions
#'
#' @param x1 An object of class \linkS4class{dph}.
#' @param x2 An object of class \linkS4class{dph}.
#' @param prob Probability for first object.
#'
#' @return An object of class \linkS4class{dph}.
#' @export
#'
#' @examples
#' dph1 <- dph(structure = "general", dimension = 3)
#' dph2 <- dph(structure = "general", dimension = 5)
#' dph_mix <- mixture(dph1, dph2, 0.5)
#' dph_mix
setMethod(
  "mixture", signature(x1 = "dph", x2 = "dph"),
  function(x1, x2, prob) {
    n1 <- length(x1@pars$alpha)
    n2 <- length(x2@pars$alpha)
    alpha <- c(prob * x1@pars$alpha, (1 - prob) * x2@pars$alpha)
    S1 <- rbind(x1@pars$S, matrix(0, n2, n1))
    S2 <- rbind(matrix(0, n1, n2), x2@pars$S)
    ph(alpha = alpha, S = cbind(S1, S2))
  }
)

#' Nfold method for phase-type distributions
#'
#' @param x1 An object of class \linkS4class{ph}.
#' @param x2 An object of class \linkS4class{dph}.
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' dph1 <- dph(structure = "general", dimension = 3)
#' dph2 <- dph(structure = "general", dimension = 2)
#' ph0 <- ph(structure = "general", dimension = 2)
#' Nfold(dph1, ph0)
#' Nfold(dph1, dph2)
setMethod(
  "Nfold", signature(x1 = "dph"),
  function(x1, x2) {
    if (methods::is(x2, "ph")) {
      exit_vect <- -x2@pars$S %*% rep(1, length(x2@pars$alpha))
      alpha <- kronecker(x1@pars$alpha, x2@pars$alpha)
      S <- kronecker(diag(length(x1@pars$alpha)), x2@pars$S) + kronecker(x1@pars$S, exit_vect %*% x2@pars$alpha)
      return(ph(alpha = alpha, S = S))
    } else if (methods::is(x2, "dph")) {
      exit_vect <- rep(1, length(x2@pars$alpha)) - x2@pars$S %*% rep(1, length(x2@pars$alpha))
      alpha <- kronecker(x1@pars$alpha, x2@pars$alpha)
      S <- kronecker(diag(length(x1@pars$alpha)), x2@pars$S) + kronecker(x1@pars$S, exit_vect %*% x2@pars$alpha)
      return(dph(alpha = alpha, S = S))
    } else {
      stop("wrong type of objects")
    }
  }
)

#' Show method for discrete phase-type distributions
#'
#' @param object An object of class \linkS4class{dph}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "dph", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@pars)
})

#' Coef method for dph Class
#'
#' @param object An object of class \linkS4class{dph}.
#'
#' @return Parameters of dph model.
#' @export
#'
#' @examples
#' obj <- dph(structure = "general", dim = 3)
#' coef(obj)
setMethod("coef", c(object = "dph"), function(object) {
  object@pars
})

#' Moment method for discrete phase-type distributions
#'
#' @param x An object of class \linkS4class{dph}.
#' @param k A positive integer (moment order).
#'
#' @return The factional moment of the \linkS4class{dph} object.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- dph(structure = "general", dimension = 3)
#' moment(obj, 2)
setMethod(
  "moment", signature(x = "dph"),
  function(x, k = 1) {
    if (k <= 0) {
      return("k should be positive")
    }
    if ((k %% 1) != 0) {
      return("k should be an integer")
    }
    m1 <- matrix_power(k - 1, x@pars$S)
    m2 <- matrix_power(k, solve(diag(nrow(x@pars$S)) - x@pars$S))
    factorial(k) * sum(x@pars$alpha %*% m1 %*% m2)
  }
)

#' Mean method for discrete phase-type distributions
#'
#' @param x An object of class \linkS4class{dph}.
#'
#' @return The raw first moment of the \linkS4class{dph} object.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- dph(structure = "general", dimension = 3)
#' mean(obj)
setMethod(
  "mean", signature(x = "dph"),
  function(x) {
    m <- solve(diag(nrow(x@pars$S)) - x@pars$S)
    sum(x@pars$alpha %*% m)
  }
)

#' Var method for discrete phase-type distributions
#'
#' @param x An object of class \linkS4class{dph}.
#'
#' @return The variance of the \linkS4class{dph} object.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- dph(structure = "general", dimension = 3)
#' var(obj)
setMethod(
  "var", signature(x = "dph"),
  function(x) {
    m <- solve(diag(nrow(x@pars$S)) - x@pars$S)
    fm <- sum(x@pars$alpha %*% m)
    m2 <- matrix_power(2, m)
    sm <- 2 * sum(x@pars$alpha %*% x@pars$S %*% m2)
    sm + fm - fm^2
  }
)

#' Pgf Method for discrete phase-type distributions
#'
#' @param x An object of class \linkS4class{dph}.
#' @param z A vector of real values.
#'
#' @return The probability generating of the \linkS4class{dph} object at the
#'  given locations.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- dph(structure = "general", dimension = 3)
#' pgf(obj, 0.5)
setMethod(
  "pgf", signature(x = "dph"),
  function(x, z) {
    if (any(abs(z) > 1)) {
      stop("z should between -1 and 1")
    }
    dph_pgf(z, x@pars$alpha, x@pars$S)
  }
)

#' Simulation method for phase-type distributions
#'
#' @param x An object of class \linkS4class{dph}.
#' @param n An integer of length of realization.
#'
#' @return A realization of independent and identically distributed discrete
#'  phase-type variables.
#' @export
#'
#' @examples
#' obj <- dph(structure = "general")
#' sim(obj, n = 100)
setMethod("sim", c(x = "dph"), function(x, n = 1000) {
  d <- length(x@pars$alpha)
  exit_vec <- 1 - rowSums(x@pars$S)
  trans_mat <- cbind(x@pars$S, exit_vec)
  aux_vec <- rep(0, d + 1)
  aux_vec[d + 1] <- 1
  trans_mat <- rbind(trans_mat, aux_vec)
  rdphasetype(n, x@pars$alpha, trans_mat)
})

#' Density method for discrete phase-type distributions
#'
#' @param x An object of class \linkS4class{dph}.
#' @param y A vector of locations.
#'
#' @return A vector containing the density evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- dph(structure = "general")
#' dens(obj, c(1, 2, 3))
setMethod("dens", c(x = "dph"), function(x, y) {
  y_inf <- (y == Inf)
  dens <- y
  dens[!y_inf] <- dphdensity(y, x@pars$alpha, x@pars$S)
  dens[y_inf] <- 0
  dens
})

#' Distribution method for discrete phase-type distributions
#'
#' @param x An object of class \linkS4class{dph}.
#' @param q A vector of locations.
#' @param lower.tail Logical parameter specifying whether lower tail (CDF) or
#' upper tail is computed.
#'
#' @return A vector containing the CDF evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- dph(structure = "general")
#' cdf(obj, c(1, 2, 3))
setMethod("cdf", c(x = "dph"), function(x,
                                        q,
                                        lower.tail = TRUE) {
  q_inf <- (q == Inf)
  cdf <- q
  cdf[!q_inf] <- dphcdf(q[!q_inf], x@pars$alpha, x@pars$S, lower.tail)
  cdf[q_inf] <- as.numeric(1 * lower.tail)
  cdf
})

#' TVR Method for dph Class
#'
#' @param x An object of class \linkS4class{dph}.
#' @param rew A vector of rewards.
#'
#' @return An object of the of class \linkS4class{dph}.
#' @export
#'
#' @examples
#' obj <- dph(structure = "general")
#' TVR(obj, c(1, 0, 1))
setMethod("TVR", c(x = "dph"), function(x, rew) {
  if (length(x@pars$alpha) != length(rew)) {
    stop("vector of rewards of wrong dimension")
  }
  if (any((rew %% 1) != 0) | any(rew < 0)) {
    stop("vector of rewards must contain only non-negative integeras")
  }
  alpha_copy <- x@pars$alpha
  S_copy <- x@pars$S

  mat_sizes <- rew + (rew == 0)

  rew_tilde <- NULL
  alpha_tilde <- NULL
  S_tilde <- NULL
  for (i in 1:length(rew)) {
    S_row <- NULL
    for (j in 1:length(rew)) {
      mat_aux <- matrix(0, nrow = mat_sizes[i], ncol = mat_sizes[j])
      if (i == j) {
        mat_aux[-mat_sizes[i], -1] <- diag(1, mat_sizes[i] - 1)
      }
      mat_aux[mat_sizes[i], 1] <- S_copy[i, j]
      S_row <- cbind(S_row, mat_aux)
    }
    S_tilde <- rbind(S_tilde, S_row)

    alpha_aux <- rep(0, mat_sizes[i])
    alpha_aux[1] <- alpha_copy[i]
    alpha_tilde <- c(alpha_tilde, alpha_aux)

    if (rew[i] != 0) {
      rew_tilde <- c(rew_tilde, rep(1, rew[i]))
    } else {
      rew_tilde <- c(rew_tilde, 0)
    }
  }

  mar_par <- tvr_dph(alpha_tilde, S_tilde, rew_tilde)
  dph(alpha = mar_par[[1]], S = mar_par[[2]])
})


#' Fit method for dph class
#'
#' @param x An object of class \linkS4class{dph}.
#' @param y Vector or data.
#' @param weight Vector of weights.
#' @param stepsEM Number of EM steps to be performed.
#' @param every Number of iterations between likelihood display updates.
#'
#' @return An object of class \linkS4class{dph}.
#'
#' @export
#'
#' @examples
#' obj <- dph(structure = "general", dimension = 2)
#' data <- sim(obj, n = 100)
#' fit(obj, data, stepsEM = 100, every = 20)
setMethod(
  "fit", c(x = "dph", y = "ANY"),
  function(x,
           y,
           weight = numeric(0),
           stepsEM = 1000,
           every = 100) {
    if (!all(y > 0)) {
      stop("data should be positive")
    }
    if (!all(weight >= 0)) {
      stop("weights should be non-negative")
    }

    A <- data_aggregation(y, weight)
    y <- A$un_obs
    weight <- A$weights

    dph_par <- x@pars
    alpha_fit <- clone_vector(dph_par$alpha)
    S_fit <- clone_matrix(dph_par$S)

    options(digits.secs = 4)
    cat(format(Sys.time(), format = "%H:%M:%OS"), ": EM started", sep = "")
    cat("\n", sep = "")

    for (k in 1:stepsEM) {
      EMstep_dph(alpha_fit, S_fit, y, weight)
      if (k %% every == 0) {
        cat("\r", "iteration:", k,
          ", logLik:", logLikelihoodDPH(alpha_fit, S_fit, y, weight),
          sep = " "
        )
      }
    }
    x@pars$alpha <- alpha_fit
    x@pars$S <- S_fit
    x@fit <- list(
      logLik = logLikelihoodDPH(alpha_fit, S_fit, y, weight),
      nobs = sum(A$weights)
    )

    cat("\n", format(Sys.time(), format = "%H:%M:%OS"), ": EM finalized", sep = "")
    cat("\n", sep = "")

    x
  }
)

#' MoE method for dph Class
#'
#' @param x An object of class \linkS4class{dph}.
#' @param formula A regression formula.
#' @param data A data frame.
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
#' x <- dph(structure = "general")
#' n <- 100
#' responses <- rpois(n, 3) + 1
#' covariate <- data.frame(age = sample(18:65, n, replace = TRUE) / 100, income = runif(n, 0, 0.99))
#' f <- responses ~ age + income # regression formula
#' MoE(x = x, formula = f, y = responses, data = covariate, stepsEM = 20)
setMethod(
  "MoE", c(x = "dph"),
  function(x,
           formula,
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
    S_fit <- clone_matrix(x@pars$S)
    c <- c()
    for (i in 1:p) c <- c(c, rep(i, n))
    extended_x <- matrix(t(as.matrix(frame[, -1])), nrow = n * p, ncol = d, byrow = TRUE)
    dm <- data.frame(Class = c, extended_x)
    names(dm)[-1] <- names(frame)[-1]
    ndm <- data.frame(dm[dm$Class == 1, -1])
    names(ndm) <- names(dm)[-1]
    for (k in 1:stepsEM) {
      B_matrix <- EMstep_dph_MoE(alpha_vecs, S_fit, frame[, 1], weight)
      wt <- reshape2::melt(B_matrix)[, 3]
      wt[wt < 1e-22] <- wt[wt < 1e-22] + 1e-22
      if (k == 1 | rand_init == TRUE) {
        multinom_model <- nnet::multinom(Class ~ ., data = dm, weights = wt, trace = F, MaxNWts = maxWts)
      } else {
        multinom_model <- nnet::multinom(Class ~ ., data = dm, weights = wt, trace = F, Wts = multinom_model$wts, MaxNWts = maxWts)
      }
      alpha_vecs <- stats::predict(multinom_model, type = "probs", newdata = ndm)
      if (k %% every == 0) {
        ll <- logLikelihoodDPH_MoE(alpha_vecs, S_fit, frame[, 1], weight)
        cat("\r", "iteration:", k, ", logLik:", ll, sep = " ")
      }
    }
    cat("\n", sep = "")
    list(alpha = alpha_vecs, S = S_fit, mm = multinom_model)
  }
)
