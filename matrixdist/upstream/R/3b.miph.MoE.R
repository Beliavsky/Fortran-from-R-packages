#' Fit method for mph/miph class, using mixture-of-experts regression
#'
#' @param x An object of class \linkS4class{mph}.
#' @param formula a regression formula.
#' @param y A matrix of observations.
#' @param data A data frame of covariates (they need to be scaled for the regression).
#' @param alpha_mat Matrix with initial distribution vectors for each row of observations.
#' @param delta Matrix with right-censoring indicators (1 uncensored, 0 right censored).
#' @param stepsEM Number of EM steps to be performed.
#' @param r Sub-sampling parameter, defaults to 1 (not supported for this method).
#' @param maxit Maximum number of iterations when optimizing the g function (inhomogeneous likelihood).
#' @param reltol Relative tolerance when optimizing g function.
#' @param rand_init Random initiation in the R-step of the EM algorithm.
#'
#' @export
#'
#' @importFrom reshape2 melt
#' @importFrom methods is
#' @importFrom stats model.frame predict optim
#' @importFrom nnet multinom
#'
#' @examples
#' under_mph <- mph(structure = c("general", "general"), dimension = 3) 
#' x <-  miph(under_mph, gfun = c("weibull", "weibull"), gfun_pars = list(c(2), c(3)))
#' n <- 100
#' responses <- cbind(rexp(n), rweibull(n, 2, 3))
#' covariates <- data.frame(age = sample(18:65, n, replace = TRUE) / 100, income = runif(n, 0, 0.99))
#' f <- responses ~ age + income
#' MoE(x = x, formula = f, y = responses, data = covariates, stepsEM = 20)
setMethod(
  "MoE", c(x = "mph", y = "ANY"),
  function(x,
           formula,
           y,
           data,
           alpha_mat = NULL,
           delta = numeric(0),
           stepsEM = 1000,
           r = 1,
           maxit = 100,
           reltol = 1e-8,
           rand_init = T) {
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

      opt_fun <- nnet_miph_LL
    }

    S_fit <- x@pars$S

    if (length(delta) == 0) {
      delta <- matrix(1, nrow(y), ncol(y))
    }

    alpha_fit <- x@pars$alpha
    p <- length(alpha_fit)
    frame <- stats::model.frame(formula, data = data)
    n <- nrow(frame)
    d2 <- ncol(frame) - 1 # number of covariates

    if (is.null(alpha_mat)) alpha_mat <- matrix(alpha_fit, ncol = p, nrow = n, byrow = TRUE) # repeats alpha n times
    c <- c()
    for (i in 1:p) c <- c(c, rep(i, n))

    extended_x <- matrix(t(as.matrix(frame[, -1])), nrow = n * p, ncol = d2, byrow = TRUE) # extended form of covariates
    dm <- data.frame(Class = c, extended_x) # data frame with all classes and covariates
    names(dm)[-1] <- names(frame)[-1]
    ndm <- data.frame(dm[dm$Class == 1, -1])
    names(ndm) <- names(dm)[-1]

    if (r < 1) {
      y_full <- y
      delta_full <- delta
      dm_full <- dm
      ndm_full <- ndm
    }

    fnn <- nnet_EM_step_mph

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
          dm <- dm_full[index, ]
          ndm <- ndm_full[index, ]
        }

        aux <- fnn(alpha_mat, S_fit, y, delta)

        B_matrix <- aux$Bmatrix
        wt <- reshape2::melt(B_matrix)[, 3]
        wt[wt < 1e-22] <- wt[wt < 1e-22] + 1e-22

        if (k == 1 | rand_init == TRUE) {
          multinom_model <- nnet::multinom(Class ~ ., data = dm, weights = wt, trace = F)
          # Class is the response and dm columns are used as predictors
        } else {
          multinom_model <- nnet::multinom(Class ~ ., data = dm, weights = wt, trace = F, Wts = multinom_model$wts)
        }
        alpha_mat <- stats::predict(multinom_model, type = "probs", newdata = ndm) # these are the new estimates for initial distribution vectors

        S_fit <- aux$S
        cat("\r", "iteration:", k,
          ", logLik:", aux$logLik,
          sep = " "
        )
      }
      x@pars$alpha <- alpha_mat # C++
      x@pars$S <- S_fit # C++

      x@fit <- list(
        # logLik = nnet_mph_LL(x, y, delta),
        loglik = sum(log(dens(x, y, delta))),
        nobs = nrow(y),
        nnet = multinom_model
      )
    }
    if (is_miph) {
      for (k in 1:stepsEM) {
        # sub-sampling
        if (r < 1) {
          index <- sample(1:nrow(y_full), size = floor(r * nrow(y_full)))

          y <- as.matrix(y_full[index, ])
          delta <- as.matrix(delta_full[index, ])
          dm <- dm_full[index, ]
          ndm <- ndm_full[index, ]
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
        aux <- fnn(alpha_mat, S_fit, trans, delta)

        S_fit <- aux$S
        B_matrix <- aux$Bmatrix
        wt <- reshape2::melt(B_matrix)[, 3]
        wt[wt < 1e-22] <- wt[wt < 1e-22] + 1e-22

        if (k == 1 | rand_init == TRUE) {
          multinom_model <- nnet::multinom(Class ~ ., data = dm, weights = wt, trace = F)
          # Class is the response and dm columns are used as predictors
        } else {
          multinom_model <- nnet::multinom(Class ~ ., data = dm, weights = wt, trace = F, Wts = multinom_model$wts)
        }
        alpha_mat <- stats::predict(multinom_model, type = "probs", newdata = ndm) # these are the new estimates for initial distribution vectors

        x@pars$alpha <- alpha_mat
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
      }
      x@pars$alpha <- alpha_mat
      x@pars$S <- S_fit
      x@gfun$pars <- par_g

      x@fit <- list(
        logLik = nnet_miph_LL(x, y, delta, x@gfun$pars),
        nobs = nrow(y),
        nnet = multinom_model
      )
    }
    cat("\n", sep = "")

    cat("\n", format(Sys.time(), format = "%H:%M:%OS"), ": EM finalized", sep = "")
    cat("\n", sep = "")

    x
  }
)

# multivariate loglikelihood to be optimized
nnet_miph_LL <- function(x,
                         obs,
                         delta,
                         gfun_pars) {
  x@gfun$pars <- gfun_pars
  res <- dens(x, y = obs, delta = delta)

  sum(log(res))
}

nnet_EM_step_mph <- function(alpha_mat, S_list, y, delta) {
  p <- ncol(alpha_mat)
  n <- nrow(y)
  d <- length(S_list)

  #---- Matrix integrals
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
  #---- Intermediate quantities
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
          a_tilde_ki[m, k, i] <- sum(alpha_mat[m, ] * a_kij[m, , i, k] * a_k_minus_i[m, , i])
        } else {
          a_tilde_ki[m, k, i] <- 0
        }
      }
    }
  }
  ###
  a <- array(NA, c(n))
  for (m in 1:n) {
    a[m] <- sum(alpha_mat[m, ] * a_k[m, ])
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
          b_ski[m, s, k, i] <- sum(alpha_mat[m, ] * a_k_minus_i[m, , i] * b_skij[m, s, k, i, ])
        }
      }
    }
  }

  #---- E and M steps
  ###
  # E STEP
  ###
  EB_k <- matrix(NA, n, p)
  for (m in 1:n) {
    for (k in 1:p) {
      EB_k[m, k] <- d * alpha_mat[m, k] * a_k[m, k] / a[m]
    }
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
  list(Bmatrix = EB_k, S = ll, logLik = sum(log(a)))
}
