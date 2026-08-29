#' Phase-type distributions
#'
#' Class of objects for phase-type distributions.
#'
#' @slot name Name of the phase-type distribution.
#' @slot pars A list comprising of the parameters.
#' @slot fit A list containing estimation information.
#'
#' @return Class object.
#' @export
#'
setClass("ph",
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

#' Constructor function for phase-type distributions
#'
#' @param alpha A probability vector.
#' @param S A sub-intensity matrix.
#' @param structure A valid ph structure: `"general"`, `"coxian"`,
#' `"hyperexponential"`, `"gcoxian"`, or `"gerlang"`.
#' @param dimension The dimension of the ph structure (if structure is provided).
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' ph(structure = "gcoxian", dimension = 5)
#' ph(alpha = c(.5, .5), S = matrix(c(-1, .5, .5, -1), 2, 2))
ph <- function(alpha = NULL, S = NULL, structure = NULL, dimension = 3) {
  if (any(is.null(alpha)) & any(is.null(S)) & is.null(structure)) {
    stop("input a vector and matrix, or a structure")
  }
  if (!is.null(structure)) {
    rs <- random_structure(dimension, structure = structure)
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
  methods::new("ph",
    name = paste(name, " ph(", length(alpha), ")", sep = ""),
    pars = list(alpha = alpha, S = S)
  )
}

#' Sum method for phase-type distributions
#'
#' @param e1 An object of class \linkS4class{ph}.
#' @param e2 An object of class \linkS4class{ph}.
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' ph1 <- ph(structure = "general", dimension = 3)
#' ph2 <- ph(structure = "gcoxian", dimension = 5)
#' ph_sum <- ph1 + ph2
#' ph_sum
setMethod(
  "+", signature(e1 = "ph", e2 = "ph"),
  function(e1, e2) {
    if (methods::is(e1, "iph") | methods::is(e2, "iph")) {
      stop("objects to be added should be ph")
    }
    L <- sum_ph(e1@pars$alpha, e1@pars$S, e2@pars$alpha, e2@pars$S)
    ph(alpha = L$alpha, S = L$S)
  }
)

kronecker_sum <- function(A, B) {
  n <- nrow(A)
  m <- nrow(B)
  kronecker(A, diag(m)) + kronecker(diag(n), B)
}

#' Minimum method for phase-type distributions
#'
#' @param x1 An object of class \linkS4class{ph}.
#' @param x2 An object of class \linkS4class{ph}.
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' ph1 <- ph(structure = "general", dimension = 3)
#' ph2 <- ph(structure = "gcoxian", dimension = 5)
#' ph_min <- minimum(ph1, ph2)
#' ph_min
setMethod(
  "minimum", signature(x1 = "ph", x2 = "ph"),
  function(x1, x2) {
    alpha <- kronecker(x1@pars$alpha, x2@pars$alpha)
    S <- kronecker_sum(x1@pars$S, x2@pars$S)
    ph(alpha = alpha, S = S)
  }
)

#' Maximum method for phase-type distributions
#'
#' @param x1 An object of class \linkS4class{ph}.
#' @param x2 An object of class \linkS4class{ph}.
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' ph1 <- ph(structure = "general", dimension = 3)
#' ph2 <- ph(structure = "gcoxian", dimension = 5)
#' ph_max <- maximum(ph1, ph2)
#' ph_max
setMethod(
  "maximum", signature(x1 = "ph", x2 = "ph"),
  function(x1, x2) {
    n1 <- length(x1@pars$alpha)
    n2 <- length(x2@pars$alpha)
    alpha <- c(kronecker(x1@pars$alpha, x2@pars$alpha), rep(0, n1 + n2))
    S1 <- rbind(kronecker_sum(x1@pars$S, x2@pars$S), matrix(0, n1 + n2, n1 * n2))
    S2 <- rbind(kronecker(diag(n1), -rowSums(x2@pars$S)), x1@pars$S, matrix(0, n2, n1))
    S3 <- rbind(kronecker(-rowSums(x1@pars$S), diag(n2)), matrix(0, n1, n2), x2@pars$S)
    ph(alpha = alpha, S = cbind(S1, S2, S3))
  }
)

#' Mixture method for phase-type distributions
#'
#' @param x1 An object of class \linkS4class{ph}.
#' @param x2 An object of class \linkS4class{ph}.
#' @param prob Probability for first object.
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' ph1 <- ph(structure = "general", dimension = 3)
#' ph2 <- ph(structure = "gcoxian", dimension = 5)
#' ph_mix <- mixture(ph1, ph2, 0.5)
#' ph_mix
setMethod(
  "mixture", signature(x1 = "ph", x2 = "ph"),
  function(x1, x2, prob) {
    n1 <- length(x1@pars$alpha)
    n2 <- length(x2@pars$alpha)
    alpha <- c(prob * x1@pars$alpha, (1 - prob) * x2@pars$alpha)
    S1 <- rbind(x1@pars$S, matrix(0, n2, n1))
    S2 <- rbind(matrix(0, n1, n2), x2@pars$S)
    ph(alpha = alpha, S = cbind(S1, S2))
  }
)

#' Moment method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#' @param k A positive integer (moment order).
#'
#' @return The raw moment of the \linkS4class{ph} (or underlying \linkS4class{ph}) object.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- ph(structure = "general", dimension = 3)
#' moment(obj, 2)
setMethod(
  "moment", signature(x = "ph"),
  function(x, k = 1) {
    if (k <= 0) {
      stop("k should be positive")
    }
    if ((k %% 1) != 0) {
      stop("k should be an integer")
    }
    if (methods::is(x, "iph")) {
      warning("moment of undelying ph structure is provided for iph objects")
    }
    prod <- matrix_power(k, solve(-x@pars$S))
    factorial(k) * sum(x@pars$alpha %*% prod)
  }
)

#' Mean method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#'
#' @return The raw first moment of the \linkS4class{ph} (or underlying \linkS4class{ph}) object.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- ph(structure = "general", dimension = 3)
#' mean(obj)
setMethod(
  "mean", signature(x = "ph"),
  function(x) {
    if (methods::is(x, "iph")) {
      warning("moment of undelying ph structure is provided for iph objects")
    }
    m <- solve(-x@pars$S)
    sum(x@pars$alpha %*% m)
  }
)

#' Var method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#'
#' @return The variance of the \linkS4class{ph} (or underlying \linkS4class{ph}) object.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- ph(structure = "general", dimension = 3)
#' var(obj)
setMethod(
  "var", signature(x = "ph"),
  function(x) {
    if (methods::is(x, "iph")) {
      warning("variance of undelying ph structure is provided for iph objects")
    }
    m <- solve(-x@pars$S)
    m2 <- matrix_power(2, m)
    2 * sum(x@pars$alpha %*% m2) - (sum(x@pars$alpha %*% m))^2
  }
)

#' Laplace method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#' @param r A vector of real values.
#'
#' @return The Laplace transform of the \linkS4class{ph}
#'  (or underlying \linkS4class{ph}) object at the given locations.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- ph(structure = "general", dimension = 3)
#' laplace(obj, 3)
setMethod(
  "laplace", signature(x = "ph"),
  function(x, r) {
    if (methods::is(x, "iph")) {
      warning("Laplace transform of undelying ph structure is provided for iph objects")
    }
    lim <- max(Re(eigen(x@pars$S)$values))
    if (any(r <= lim)) {
      stop("r should be above the largest real eigenvalue of S")
    }
    ph_laplace(r, x@pars$alpha, x@pars$S)
  }
)

#' Mgf method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#' @param r A vector of real values.
#'
#' @return The mgf of the \linkS4class{ph} (or underlying \linkS4class{ph}) object
#'  at the given locations.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- ph(structure = "general", dimension = 3)
#' mgf(obj, 0.4)
setMethod(
  "mgf", signature(x = "ph"),
  function(x, r) {
    if (methods::is(x, "iph")) {
      warning("mgf of undelying ph structure is provided for iph objects")
    }
    lim <- -max(Re(eigen(x@pars$S)$values))
    if (any(r > lim)) {
      stop("r should be below the negative largest real eigenvalue of S")
    }
    ph_laplace(-r, x@pars$alpha, x@pars$S)
  }
)

#' Show method for phase-type distributions
#'
#' @param object An object of class \linkS4class{ph}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "ph", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@pars)
})

#' Simulation method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#' @param n An integer of length of realization.
#'
#' @return A realization of independent and identically distributed phase-type variables.
#' @export
#'
#' @examples
#' obj <- ph(structure = "general")
#' sim(obj, n = 100)
setMethod("sim", c(x = "ph"), function(x, n = 1000) {
  rphasetype(n, x@pars$alpha, x@pars$S)
})

#' Density method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#' @param y A vector of locations.
#'
#' @return A vector containing the density evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- ph(structure = "general")
#' dens(obj, c(1, 2, 3))
setMethod("dens", c(x = "ph"), function(x, y) {
  y_inf <- (y == Inf)
  dens <- y
  dens[!y_inf] <- phdensity(y, x@pars$alpha, x@pars$S)
  dens[y_inf] <- 0
  dens
})

#' Distribution method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#' @param q A vector of locations.
#' @param lower.tail Logical parameter specifying whether lower tail (CDF) or
#' upper tail is computed.
#'
#' @return A vector containing the CDF evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- ph(structure = "general")
#' cdf(obj, c(1, 2, 3))
setMethod("cdf", c(x = "ph"), function(x,
                                       q,
                                       lower.tail = TRUE) {
  q_inf <- (q == Inf)
  cdf <- q
  cdf[!q_inf] <- phcdf(q[!q_inf], x@pars$alpha, x@pars$S, lower.tail)
  cdf[q_inf] <- as.numeric(1 * lower.tail)
  cdf
})

#' Hazard rate method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#' @param y A vector of locations.
#'
#' @return A vector containing the hazard rate evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- ph(structure = "general")
#' haz(obj, c(1, 2, 3))
setMethod("haz", c(x = "ph"), function(x, y) {
  d <- dens(x, y)
  s <- cdf(x, y, lower.tail = FALSE)
  d / s
})

#' Quantile method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#' @param p A vector of probabilities.
#'
#' @return A vector containing the quantile evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- ph(structure = "general")
#' quan(obj, c(0.5, 0.9, 0.99))
setMethod("quan", c(x = "ph"), function(x,
                                        p) {
  quan <- numeric(length(p))
  for (i in seq_along(p)) {
    quan[i] <- stats::uniroot(f = function(q) p[i] - cdf(x, 1 / (1 - q) - 1), interval = c(0, 1))$root
  }
  1 / (1 - quan) - 1
})

#' Fit method for ph class
#'
#' @param x An object of class \linkS4class{ph}.
#' @param y Vector or data.
#' @param weight Vector of weights.
#' @param rcen Vector of right-censored observations.
#' @param rcenweight Vector of weights for right-censored observations.
#' @param stepsEM Number of EM steps to be performed.
#' @param methods Methods to use for matrix exponential calculation: RM, UNI or PADE.
#' @param rkstep Runge-Kutta step size (optional).
#' @param uni_epsilon Epsilon parameter for uniformization method.
#' @param maxit Maximum number of iterations when optimizing g function.
#' @param reltol Relative tolerance when optimizing g function.
#' @param every Number of iterations between likelihood display updates.
#' @param r Sub-sampling proportion for stochastic EM, defaults to 1.
#'
#' @return An object of class \linkS4class{ph}.
#'
#' @importFrom grDevices dev.off
#' @importFrom graphics hist legend lines
#' @importFrom utils head
#'
#' @export
#'
#' @examples
#' obj <- iph(ph(structure = "general", dimension = 2), gfun = "weibull", gfun_pars = 2)
#' data <- sim(obj, n = 100)
#' fit(obj, data, stepsEM = 100, every = 20)
setMethod(
  "fit", c(x = "ph", y = "ANY"),
  function(x,
           y,
           weight = numeric(0),
           rcen = numeric(0),
           rcenweight = numeric(0),
           stepsEM = 1000,
           methods = c("RK", "RK"),
           rkstep = NA,
           uni_epsilon = NA,
           maxit = 100,
           reltol = 1e-8,
           every = 100,
           r = 1) {
    EMstep <- eval(parse(text = paste("EMstep_", methods[1], sep = "")))
    if (!all(c(y, rcen) > 0)) {
      stop("data should be positive")
    }
    if (!all(c(weight, rcenweight) >= 0)) {
      stop("weights should be non-negative")
    }
    if (r < 1 && any(methods == "RK")) {
      stop("sub-sampling available for UNI and PADE methods")
    }
    is_iph <- methods::is(x, "iph")
    if (is_iph) {
      par_g <- x@gfun$pars
      inv_g <- x@gfun$inverse
    }
    LL <- if (is_iph) {
      eval(parse(text = paste("logLikelihoodM", x@gfun$name, "_", methods[2], sep = "")))
    } else {
      eval(parse(text = paste("logLikelihoodPH_", methods[2], sep = "")))
    }
    A <- data_aggregation(y, weight)
    y <- A$un_obs
    weight <- A$weights
    if (length(rcen) > 0) {
      B <- data_aggregation(rcen, rcenweight)
      rcen <- B$un_obs
      rcenweight <- B$weights
    }

    ph_par <- x@pars
    alpha_fit <- clone_vector(ph_par$alpha)
    S_fit <- clone_matrix(ph_par$S)

    if (r < 1) {
      y_full <- y
      weight_full <- weight
      rcen_full <- rcen
      rcenweight_full <- rcenweight
    }

    options(digits.secs = 4)
    cat(format(Sys.time(), format = "%H:%M:%OS"), ": EM started", sep = "")
    cat("\n", sep = "")

    if (!is_iph) {
      for (k in 1:stepsEM) {
        if (r < 1) {
          indices <- sample(1:length(y_full), size = floor(r * length(y_full)))
          y <- y_full[indices]
          weight <- weight_full[indices]
          if (length(rcen_full) > 0) {
            rcen <- rcen_full[indices]
            rcenweight <- rcenweight_full[indices]
          }
        }

        epsilon1 <- switch(which(methods[1] == c("RK", "UNI", "PADE")),
          if (!is.na(rkstep)) rkstep else default_step_length(S_fit),
          if (!is.na(uni_epsilon)) uni_epsilon else 1e-4,
          0
        )
        epsilon2 <- switch(which(methods[2] == c("RK", "UNI", "PADE")),
          if (!is.na(rkstep)) rkstep else default_step_length(S_fit),
          if (!is.na(uni_epsilon)) uni_epsilon else 1e-4,
          0
        )
        EMstep(epsilon1, alpha_fit, S_fit, y, weight, rcen, rcenweight)
        if (k %% every == 0) {
          cat("\r", "iteration:", k,
            ", logLik:", LL(epsilon2, alpha_fit, S_fit, y, weight, rcen, rcenweight),
            sep = " "
          )
        }
      }
      x@pars$alpha <- alpha_fit
      x@pars$S <- S_fit
      x@fit <- list(
        logLik = LL(epsilon2, alpha_fit, S_fit, y, weight, rcen, rcenweight),
        nobs = sum(A$weights)
      )
    }
    if (is_iph) {
      trans_weight <- weight
      trans_rcenweight <- rcenweight
      for (k in 1:stepsEM) {
        if (r < 1) {
          indices <- sample(1:length(y_full), size = floor(r * length(y_full)))
          y <- y_full[indices]
          weight <- weight_full[indices]
          if (length(rcen_full) > 0) {
            rcen <- rcen_full[indices]
            rcenweight <- rcenweight_full[indices]
          }
        }

        if (x@gfun$name != "gev") {
          trans <- inv_g(par_g, y)
          trans_cens <- inv_g(par_g, rcen)
        } else {
          t <- inv_g(par_g, y, weight)
          tc <- inv_g(par_g, rcen, rcenweight)
          trans <- t$obs
          trans_weight <- t$weight
          trans_cens <- tc$obs
          trans_rcenweight <- tc$weight
        }
        epsilon1 <- switch(which(methods[1] == c("RK", "UNI", "PADE")),
          if (!is.na(rkstep)) rkstep else default_step_length(S_fit),
          if (!is.na(uni_epsilon)) uni_epsilon else 1e-4,
          0
        )
        epsilon2 <- switch(which(methods[2] == c("RK", "UNI", "PADE")),
          if (!is.na(rkstep)) rkstep else default_step_length(S_fit),
          if (!is.na(uni_epsilon)) uni_epsilon else 1e-4,
          0
        )
        EMstep(epsilon1, alpha_fit, S_fit, trans, trans_weight, trans_cens, trans_rcenweight)
        opt <- suppressWarnings(
          stats::optim(
            par = par_g,
            fn = LL,
            h = epsilon2,
            alpha = alpha_fit,
            S = S_fit,
            obs = y,
            weight = weight,
            rcens = rcen,
            rcweight = rcenweight,
            hessian = FALSE,
            control = list(
              maxit = maxit,
              reltol = reltol,
              fnscale = -1
            )
          )
        )
        par_g <- opt$par
        if (k %% every == 0) {
          cat("\r", "iteration:", k,
            ", logLik:", opt$value,
            sep = " "
          )
        }
      }
      x@pars$alpha <- alpha_fit
      x@pars$S <- S_fit
      x@fit <- list(
        logLik = opt$value,
        nobs = sum(A$weights)
      )
      x <- iph(x, gfun = x@gfun$name, gfun_pars = par_g)
    }

    cat("\n", format(Sys.time(), format = "%H:%M:%OS"), ": EM finalized", sep = "")
    cat("\n", sep = "")

    x
  }
)

data_aggregation <- function(y, w) {
  if (length(w) == 0) w <- rep(1, length(y))
  observations <- cbind(y, w)
  mat <- data.frame(observations)
  names(mat) <- c("obs", "weight")
  if (length(y) == 0) {
    return(list(un_obs = numeric(0), weights = numeric(0)))
  } else {
    agg <- stats::aggregate(mat$weight,
      by = list(un_obs = mat$obs),
      FUN = sum
    )
    return(list(un_obs = agg$un_obs, weights = agg$x))
  }
}

#' Loglikelihood method for ph class
#'
#' @param object An object of class \linkS4class{ph}.
#'
#' @return An object of class logLik.
#' @export
#'
#' @examples
#' obj <- iph(ph(structure = "general", dimension = 2), gfun = "weibull", gfun_pars = 2)
#' data <- sim(obj, n = 100)
#' fitted_ph <- fit(obj, data, stepsEM = 10)
#' logLik(fitted_ph)
setMethod("logLik", "ph", function(object) {
  ll <- object@fit$logLik
  attr(ll, "nobs") <- object@fit$nobs
  attr(ll, "df") <- sum(unlist(coef(object)) != 0) - 1
  class(ll) <- "logLik"
  ll
})

#' Coef method for ph class
#'
#' @param object An object of class \linkS4class{ph}.
#'
#' @return Parameters of ph model.
#' @export
#'
#' @examples
#' obj <- ph(structure = "general")
#' coef(obj)
setMethod("coef", c(object = "ph"), function(object) {
  object@pars
})

#' LRT method for ph class
#'
#' @param x,y Objects of class \linkS4class{ph}.
#'
#' @return LRT between the models.
#' @export
#' @importFrom stats pchisq
#'
setMethod("LRT", c(x = "ph", y = "ph"), function(x, y) {
  LR <- 2 * abs(logLik(y) - logLik(x))
  degrees <- abs(attributes(logLik(y))$df - attributes(logLik(x))$df)
  c(LR = LR, p.val = pchisq(LR, df = degrees, lower.tail = FALSE))
})

#' TVR method for ph class
#'
#' @param x An object of class \linkS4class{ph}.
#' @param rew A vector of rewards.
#'
#' @return An object of the of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' obj <- ph(structure = "general")
#' TVR(obj, c(1, 2, 3))
setMethod("TVR", c(x = "ph"), function(x, rew) {
  if (length(x@pars$alpha) != length(rew)) {
    stop("vector of rewards of wrong dimension")
  }
  if (any(rew < 0)) {
    stop("vector of rewards with negative entries")
  }
  mar_par <- tvr_ph(x@pars$alpha, x@pars$S, rew)
  ph(alpha = mar_par[[1]], S = mar_par[[2]])
})
