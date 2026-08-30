
# mvnorm

### allow more than one distribution
mvnorm <- function(mean, invcholmean, chol, invchol) {

  # mvnorm chol invchol
  
  if (missing(chol) && missing(invchol))
      chol <- as.chol(ltMatrices(1, diag = TRUE))
  stopifnot(xor(missing(chol), missing(invchol)))

  if (!missing(chol)) {
      if (!is.ltMatrices(chol))
          chol <- as.ltMatrices(chol)
      scale <- as.chol(chol)
  }

  if (!missing(invchol)) {
      if (!is.ltMatrices(invchol))
          invchol <- as.ltMatrices(invchol)
      scale <- as.invchol(invchol)
  }
  ret <- list(scale = scale)
  
  if (!missing(invcholmean)) {
      stopifnot(missing(mean))
      ### note: this is not really necessary once all ls functions
      ### are able to deal with invcholmean
      # mvnorm invcholmean
      
      if (!missing(invcholmean)) {
          stopifnot(is.numeric(invcholmean))
          stopifnot(NROW(invcholmean) == dim(scale)[2L])
          if (!is.matrix(invcholmean)) {
              invcholmean <- matrix(invcholmean, nrow = NROW(invcholmean))
              rownames(invcholmean) <- names(invcholmean)
          }
          nm <- dimnames(scale)[[2L]]
          if (is.null(rownames(invcholmean)))
              rownames(invcholmean) <- nm
          if (!isTRUE(all.equal(rownames(invcholmean), nm)))
              stop("rownames of invcholmean do not match")        
          nm <- dimnames(scale)[[1L]]
          if (!is.null(nm) && dim(scale)[[2L]] == ncol(invcholmean)) {
              if (is.null(colnames(invcholmean)))
                  colnames(invcholmean) <- nm
              if (!isTRUE(all.equal(colnames(invcholmean), nm)))
                  stop("colnames of invcholmean do not match")        
          }
          ret$invcholmean <- invcholmean
      }
      
  } else {
      # mvnorm mean
      
      if (!missing(mean)) {
          stopifnot(is.numeric(mean))
          stopifnot(NROW(mean) == dim(scale)[2L])
          if (!is.matrix(mean)) {
              mean <- matrix(mean, nrow = NROW(mean))
              rownames(mean) <- names(mean)
          }
          nm <- dimnames(scale)[[2L]]
          if (is.null(rownames(mean)))
              rownames(mean) <- nm
          if (!isTRUE(all.equal(rownames(mean), nm)))
              stop("rownames of mean do not match")        
          nm <- dimnames(scale)[[1L]]
          if (!is.null(nm) && dim(scale)[[2L]] == ncol(mean)) {
              if (is.null(colnames(mean)))
                  colnames(mean) <- nm
              if (!isTRUE(all.equal(colnames(mean), nm)))
                  stop("colnames of mean do not match")        
          }
          ret$mean <- mean
      }
      
  }
  class(ret) <- "mvnorm"
  return(ret)
}

# mean2invcholmean

.m2i <- function(object) {
    if (is.chol(object$scale))
        return(solve(object$scale, object$mean))
    return(object$scale %*% object$mean)
}

# invcholmean2mean

.i2m <- function(object) {
    if (is.chol(object$scale))
        return(object$scale %*% object$invcholmean)
    return(solve(object$scale, object$invcholmean))
}

# mvnorm methods

names.mvnorm <- function(x)
    dimnames(x$scale)[[2L]]

aperm.mvnorm <- function(a, perm, ...) {

    if (!is.null(a$invcholmean))
        a$mean <- .i2m(a)
    ret <- list(scale = aperm(a$scale, perm = perm, ...))
    if (!is.null(a$mean))
        ret$mean <- a$mean[perm,,drop = FALSE]
    if (!is.null(a$invcholmean)) {
        ret$invcholmean <- .m2i(ret)
        ret$mean <- NULL
    }
    class(ret) <- "mvnorm"
    ret
}

# mvnorm simulate

simulate.mvnorm <- function(object, nsim = dim(object$scale)[1L], seed = NULL, 
                            standardize = FALSE, as.data.frame = FALSE, ...) {

    # init random seed, reset on exit
    
    ### from stats:::simulate.lm
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    

    J <- dim(object$scale)[2L]
    N <- dim(object$scale)[1L]
    if (N > 1)
        stopifnot(nsim == N)
    if (standardize) {
        if (is.chol(object$scale)) {
            object$scale <- standardize(chol = object$scale)
        } else {
            object$scale <- standardize(invchol = object$scale)
        }
    }
    Z <- matrix(rnorm(nsim * J), nrow = J)
    if (!is.null(object$invcholmean))
        Z <- Z + c(object$invcholmean)
    if (is.chol(object$scale)) {
        Y <- Mult(object$scale, Z)
    } else {
        Y <- solve(object$scale, Z)
    }
    ret <- Y
    if (!is.null(object$mean))
        ret <- ret + c(object$mean)
    rownames(ret) <- dimnames(object$scale)[[2L]]
    if (!as.data.frame)
        return(ret)
    return(as.data.frame(t(ret)))
}

# mvnorm margDist

margDist <- function(object, which, ...)
    UseMethod("margDist")

margDist.mvnorm <- function(object, which, ...) {

    if (is.chol(object$scale)) {
        ret <- list(scale = as.chol(marg_mvnorm(chol = object$scale, 
                                                which = which)$chol))
    } else {
        ret <- list(scale = as.invchol(marg_mvnorm(invchol = object$scale, 
                                                   which = which)$invchol))
    }

    if (!is.null(object$invcholmean))
        object$mean <- .i2m(object)
    if (!is.null(object$mean))
        ret$mean <- object$mean[which,,drop = FALSE]
    if (!is.null(object$invcholmean)) {
        ret$invcholmean <- .m2i(ret)
        ret$mean <- NULL
    }

    class(ret) <- "mvnorm"
    return(ret)
}

# mvnorm condDist

condDist <- function(object, which_given, given, ...)
    UseMethod("condDist")

condDist.mvnorm <- function(object, which_given = 1L, given, ...) {

    if (!is.null(object$invcholmean))
        object$mean <- .i2m(object)
    if (!is.null(object$mean)) {
        if (ncol(object$mean) == 1L && NCOL(given) > 1L) {
            given <- given - object$mean[which_given,,drop = TRUE]
        } else {
            given <- c(given) - object$mean[which_given,,drop = FALSE]
        }
    }

    if (is.chol(object$scale)) {
        ret <- cond_mvnorm(chol = object$scale, which_given = which_given, 
                           given = given, ...)
        ret$scale <- as.chol(ret$chol)
        ret$chol <- NULL
    } else {
        ret <- cond_mvnorm(invchol = object$scale, which_given = which_given, 
                           given = given, ...)
        ret$scale <- as.invchol(ret$invchol)
        ret$invchol <- NULL
    }
    if (!is.null(object$mean)) {
        if (is.character(which_given)) 
            which_given <- match(which_given, dimnames(object$scale)[[2L]])
        if (ncol(object$mean) > 1L && ncol(ret$mean) > 1) {
            ret$mean <- object$mean[-which_given,,drop = FALSE] + ret$mean
        } else if (ncol(object$mean) == 1L && ncol(ret$mean) > 1L) {
            ret$mean <- object$mean[-which_given,,drop = TRUE] + ret$mean
        } else {
            ret$mean <- object$mean[-which_given,,drop = FALSE] + c(ret$mean)
        }
    }
    if (!is.null(object$invcholmean)) {
        ret$invcholmean <- .m2i(ret)
        ret$mean <- NULL
    }
    class(ret) <- "mvnorm"
    return(ret)
}

# mvnorm logLik

logLik.mvnorm <- function(object, obs, lower, upper, standardize = FALSE, 
                          ...) {
    # argchecks
    
    args <- c(object, list(...))
    nargs <- missing(obs) + missing(lower) + missing(upper)
    stopifnot(nargs < 3L)

    nmobs <- NULL
    if (!missing(obs)) {
        if (!is.null(obs)) {
            stopifnot(is.matrix(obs))
            stopifnot(!is.null(rownames(obs)))
            nmobs <- rownames(obs)
        }
    }
    nmlower <- nmupper <- nmlu <- NULL
    if (!missing(lower)) {
        if (!is.null(lower)) {
            stopifnot(is.matrix(lower))
            stopifnot(!is.null(rownames(lower)))
            nmlu <- nmlower <- rownames(lower)
        }
    }
    if (!missing(upper)) {
        if (!is.null(lower)) {
            stopifnot(is.matrix(upper))
            stopifnot(!is.null(rownames(upper)))
            nmupper <- rownames(upper)
            if (!missing(lower)) {
                stopifnot(isTRUE(all.equal(nmlower, nmupper)))
            } else {
                nmlu <- nmupper
            }
        }
    }

    nm <- c(nmobs, nmlu)
    no <- names(object)
    ### we allow dimensions w/o data
    stopifnot(all(nm %in% no))
    perm <- NULL
    if (!isTRUE(all.equal(nm, no)))
        perm <- c(nm, no[!no %in% nm])

    if (!missing(obs)) args$obs <- obs
    if (!missing(lower)) args$lower <- lower
    if (!missing(upper)) args$upper <- upper
    

    # logLik standardize
    
    if (is.chol(object$scale)) {
        names(args)[names(args) == "scale"] <- "chol"
        if (standardize)
           args$chol <- object$scale <- standardize(chol = args$chol)
    } else {
        names(args)[names(args) == "scale"] <- "invchol"
        if (standardize)
            args$invchol <- object$scale <- standardize(invchol = args$invchol)
    }
    

    if (!is.null(perm)) {
        ll <- if (!is.null(args$logLik)) args$logLik else TRUE
        ret <- 0
        ### integrate out dimensions w/o data
        if (length(nm) < length(no) && !is.null(nmlu))
            object <- margDist(object, which = nm)
        ### continuous
        if (!is.null(nmobs)) {
            ret <- ret + logLik(margDist(object, which = nmobs), 
                                obs = obs, logLik = ll) 
            ### set-up conditional distribution if there are any intervals
            if (!is.null(nmlu))
                object <- condDist(object, which_given = nmobs, given = obs)
        }
        ### interval given continuous
        if (!is.null(nmlu))
            ret <- ret + logLik(object, lower = lower, upper = upper, ...)
        return(ret)
    }

    ### base likelihood on mean for mixed data
    if (!is.null(args$invcholmean) &&
        (!missing(obs)) + (!missing(lower) || !missing(upper)) == 2L  ### mixed data
        ) {
        args$mean <- .i2m(object)
        args$invcholmean <- NULL
    }
    return(do.call("ldpmvnorm", args))
}

# mvnorm lLgrad

lLgrad <- function(object, ...)
    UseMethod("lLgrad")

lLgrad.mvnorm <- function(object, obs, lower, upper, standardize = FALSE, 
                          ...) {
    # argchecks
    
    args <- c(object, list(...))
    nargs <- missing(obs) + missing(lower) + missing(upper)
    stopifnot(nargs < 3L)

    nmobs <- NULL
    if (!missing(obs)) {
        if (!is.null(obs)) {
            stopifnot(is.matrix(obs))
            stopifnot(!is.null(rownames(obs)))
            nmobs <- rownames(obs)
        }
    }
    nmlower <- nmupper <- nmlu <- NULL
    if (!missing(lower)) {
        if (!is.null(lower)) {
            stopifnot(is.matrix(lower))
            stopifnot(!is.null(rownames(lower)))
            nmlu <- nmlower <- rownames(lower)
        }
    }
    if (!missing(upper)) {
        if (!is.null(lower)) {
            stopifnot(is.matrix(upper))
            stopifnot(!is.null(rownames(upper)))
            nmupper <- rownames(upper)
            if (!missing(lower)) {
                stopifnot(isTRUE(all.equal(nmlower, nmupper)))
            } else {
                nmlu <- nmupper
            }
        }
    }

    nm <- c(nmobs, nmlu)
    no <- names(object)
    ### we allow dimensions w/o data
    stopifnot(all(nm %in% no))
    perm <- NULL
    if (!isTRUE(all.equal(nm, no)))
        perm <- c(nm, no[!no %in% nm])

    if (!missing(obs)) args$obs <- obs
    if (!missing(lower)) args$lower <- lower
    if (!missing(upper)) args$upper <- upper
    

    ### base likelihood on mean when necessary
    ### note that post processing the mean score is time consuming
    if (!is.null(args$invcholmean) &&
        (!is.null(perm) ||      ### permutations        
         (!missing(obs)) + (!missing(lower) || !missing(upper)) == 2L  ### mixed data
        )) {
        args$mean <- .i2m(object)
        args$invcholmean <- NULL
    }

    if (is.chol(object$scale)) {
        # lLgrad chol
        
        names(args)[names(args) == "scale"] <- "chol"
        sc <- args$chol
        if (standardize)
            args$chol <- sc <- standardize(chol = args$chol)
        if (!is.null(perm)) {
            if (!attr(args$chol, "diag")) {
                diagonals(args$chol) <- 1
                sc <- args$chol
            }
            args$chol <- pc <- aperm(as.chol(args$chol), perm = perm)
            if (length(nm) < length(no))
                args$chol <- marg_mvnorm(chol = args$chol, which = nm)$chol
            args$mean <- args$mean[nm,,drop = FALSE]
        }
        ret <- do.call("sldpmvnorm", args)
        # lLgrad mean

        ### sldmvnorm returns mean score as -obs
        if (is.null(ret$mean)) ret$mean <- - ret$obs
        
        # lLgrad marginalisation

        om <- length(no) - length(nm)
        if (om > 0) {
            am <- matrix(0, nrow = om, ncol = ncol(ret$mean))
            rownames(am) <- no[!no %in% nm]
            ret$mean <- rbind(ret$mean, am)
            Jo <- dim(object$scale)[[2L]]
            pJ <- dim(args$invchol)[[2L]]
            am <- matrix(0, nrow = Jo * (Jo + 1) / 2 - pJ * (pJ + 1) / 2, 
                         ncol = dim(ret$invchol)[1L])
            byrow_orig <- attr(ret$chol, "byrow")
            ret$chol <- ltMatrices(ret$chol, byrow = TRUE)
            ### rbind only works for byrow = TRUE
            ret$chol <- ltMatrices(rbind(unclass(ret$chol), am), 
                                   byrow = TRUE, 
                                   diag = TRUE,
                                   names = perm)
            ret$chol <- ltMatrices(ret$chol, byrow = byrow_orig)
        }
        
        # lLgrad deperma

        if (!is.null(perm))
            ret$chol <- deperma(chol = sc, permuted_chol = pc, 
                                perm = match(perm, no), 
                                score_schol = ret$chol)
        
        # lLgrad destandarized

        if (standardize)
            ret$chol <- destandardize(chol = object$scale, 
                                      score_schol = ret$chol)
        
        # lLgrad diagonals

        if (!attr(sc, "diag"))
            ret$chol <- ltMatrices(Lower_tri(ret$chol, diag = FALSE),
                                   diag = FALSE, 
                                   byrow = attr(ret$chol, "byrow"), 
                                   names = dimnames(ret$chol)[[2L]])
        
        # lLgrad return

        ret$scale <- ret$chol
        ret$chol <- NULL
        if (!is.null(ret$mean)) ret$mean <- ret$mean[no,, drop = FALSE]
        ### sldpmvnorm for mixed data returns mean score only
        if (!is.null(object$invcholmean) && is.null(ret$invcholmean)) {
            J <- dim(sc)[2L]
            M <- matrix(seq_len(J^2), nrow = J, byrow = FALSE)
            idx <- M[.lt(J, diag = TRUE)]

            X <- ret$mean
            Y <- matrix(object$invcholmean, nrow = nrow(X), ncol = ncol(X))
            A <- X[rep(1:nrow(X), times = nrow(X)),,drop = FALSE] * 
                 Y[rep(1:nrow(Y), each = nrow(X)),,drop = FALSE]

            scC <- ltMatrices(A[idx,,drop = FALSE], diag = TRUE, byrow = FALSE)
            if (!attr(sc, "diag"))
                scC <- ltMatrices(Lower_tri(scC, diag = FALSE), diag = FALSE, byrow = FALSE)
            scC <- ltMatrices(scC, byrow = attr(sc, "byrow"))
            ret$scale <- ltMatrices(unclass(scC) + unclass(ret$scale), 
                                    diag = attr(sc, "diag"),
                                    byrow = attr(sc, "byrow"),
                                    names = dimnames(sc)[[2L]])
            ret$invcholmean <- Mult(sc, X, transpose = TRUE)
            ret$mean <- NULL
        }
        return(ret)
        
        
    }
    # lLgrad invchol
    
    names(args)[names(args) == "scale"] <- "invchol"
    si <- args$invchol
    if (standardize)
        args$invchol <- si <- standardize(invchol = args$invchol)
    if (!is.null(perm)) {
        if (!attr(args$invchol, "diag")) {
            diagonals(args$invchol) <- 1
            si <- args$invchol
        }
        args$invchol <- pi <- aperm(as.invchol(args$invchol), perm = perm)
        if (length(nm) < length(no))
            args$invchol <- marg_mvnorm(invchol = args$invchol,
                                        which = nm)$invchol
        args$mean <- args$mean[nm,,drop = FALSE]
    }
    ret <- do.call("sldpmvnorm", args)
    ### sldmvnorm returns mean score as -obs
    ### return only if object had mean specified
    if (is.null(ret$mean) && !is.null(args$mean)) 
        ret$mean <- - ret$obs

    # lLgrad invchol marginalisation

    om <- length(no) - length(nm)
    if (om > 0) {
        am <- matrix(0, nrow = om, ncol = ncol(ret$mean))
        rownames(am) <- no[!no %in% nm]
        ret$mean <- rbind(ret$mean, am)
        Jo <- dim(object$scale)[[2L]]
        pJ <- dim(args$invchol)[[2L]]
        am <- matrix(0, nrow = Jo * (Jo + 1) / 2 - pJ * (pJ + 1) / 2, 
                     ncol = dim(ret$invchol)[1L])
        byrow_orig <- attr(ret$invchol, "byrow")
        ret$invchol <- ltMatrices(ret$invchol, byrow = TRUE)
        ### rbind only works for byrow = TRUE
        ret$invchol <- ltMatrices(rbind(unclass(ret$invchol), am), 
                                  byrow = TRUE,
                                  diag = TRUE,
                                  names = perm)
        ret$invchol <- ltMatrices(ret$invchol, byrow = byrow_orig)
    }
    

    if (!is.null(perm))
        ret$invchol <- deperma(invchol = si, permuted_invchol = pi, 
                               perm = match(perm, no), 
                               score_schol = -vectrick(pi, ret$invchol))
    if (standardize)
        ret$invchol <- destandardize(invchol = object$scale, 
                                     score_schol = -vectrick(si, ret$invchol))
    if (!attr(si, "diag"))
        ret$invchol <- ltMatrices(Lower_tri(ret$invchol, diag = FALSE),
                                  diag = FALSE, 
                                  byrow = attr(ret$invchol, "byrow"), 
                                  names = dimnames(ret$invchol)[[2L]])
    ret$scale <- ret$invchol
    ret$invchol <- NULL
    if (!is.null(ret$mean)) ret$mean <- ret$mean[no,, drop = FALSE]
    # lLgrad invchol post

    ### sldpmvnorm for mixed data returns mean score only
    if (!is.null(object$invcholmean) && is.null(ret$invcholmean)) {
        J <- dim(si)[2L]
        M <- matrix(seq_len(J^2), nrow = J, byrow = FALSE)
        idx <- M[.lt(J, diag = TRUE)]

        X <- ret$mean
        Y <- matrix(object$invcholmean, nrow = nrow(X), ncol = ncol(X))
        A <- X[rep(1:nrow(X), times = nrow(X)),,drop = FALSE] * 
             Y[rep(1:nrow(Y), each = nrow(X)),,drop = FALSE]

        scL <- - vectrick(solve(si), A)
        scL <- ltMatrices(scL[idx,,drop = FALSE], diag = TRUE, byrow = FALSE)
        if (!attr(si, "diag"))
            scL <- ltMatrices(Lower_tri(scL, diag = FALSE), diag = FALSE, byrow = FALSE)
        scL <- ltMatrices(scL, byrow = attr(si, "byrow"))
        ret$scale <- ltMatrices(unclass(scL) + unclass(ret$scale), 
                                  diag = attr(si, "diag"),
                                  byrow = attr(si, "byrow"),
                                  names = dimnames(si)[[2L]])
        ret$invcholmean <- solve(si, X, transpose = TRUE)
        ret$mean <- FALSE
    }
    
    return(ret)
    
}

# mvnorm marginal means

mean.mvnorm <- function(x, ...) {
    if (!is.null(x$mean)) return(x$mean)
    if (!is.null(x$invcholmean)) return(.i2m(x))
    return(0)
}

# mvnorm vcov

vcov.mvnorm <- function(object, ...) {
    if (is.invchol(object$scale)) return(invchol2cov(object$scale))
    return(chol2cov(object$scale))
}

# mvnorm regression coefs

coef.mvnorm <- function(object, which = dim(object$scale)[2L], ...) {

    if (dim(object$scale)[[2L]] == 1) {
        if (!is.invchol(object$scale)) object$scale <- chol2invchol(object$scale)
        if (!is.null(object$invcholmean)) object$mean <- mean(object)
        if (is.null(object$mean)) return(NULL)
        ret <- object$mean
        attr(ret, "sigma") <- 1 / c(unclass(object$scale))
        return(ret)
    }

    if (is.character(which)) 
        which <- which(dimnames(object$scale)[[2L]] == which)

    stopifnot((length(which) == 1L) && 
              (which > 1) && 
              (which <= dim(object$scale)[2L]))

    if (!is.invchol(object$scale)) object$scale <- chol2invchol(object$scale)

    diag <- attr(object$scale, "diag")
    j <- which - 1L
    start <- j * (j + 1) / 2
    idx <- start + seq_len(which)
    ret <- Lower_tri(ltMatrices(object$scale, byrow = TRUE, diag = diag), 
                     diag = diag)[idx,,drop = FALSE]
    rownames(ret) <- dimnames(object$scale)[[2L]][seq_len(which)]
    if (!is.null(nm <- dimnames(object$scale)[[1L]]))
        colnames(ret) <- nm
    dg <- ret[nrow(ret),]
    cf <- - ret[-nrow(ret),,drop = FALSE]
    if (!is.null(object$mean)) object$invcholmean <- .m2i(object)
    if (!is.null(object$invcholmean))
        cf <- rbind("(Intercept)" = object$invcholmean[which,,drop = TRUE], cf)
    cf <- cf / matrix(dg, nrow = nrow(cf), ncol = ncol(cf), byrow = TRUE)
    cf <- cf[,,drop = TRUE]
    attr(cf, "sigma") <- 1 / dg
    return(cf)
}

