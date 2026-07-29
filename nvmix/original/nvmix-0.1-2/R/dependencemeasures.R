##' @title Dependence measures for grouped nvmix copulas 
##' 
##' @param scale n-vector giving the 'scale' parameters
##' @param qmix see ?pgnvmix()
##' @param method string ("kendall" or "spearman") which measure is computed 
##' @param groupings either rep(1, 2) (ungrouped) or c(1, 2) (grouped)
##' @param ellip.kendall logical if formula for Kendall's tau for elliptical
##'        copulas shall be used 
##' @param control list(), see get_set_param()
##' @param verbose logical (or integer: 0 = FALSE, 1 = TRUE, 2 = more output)
##'        indicating whether a warning is given if the required precision
##'        'control$pnvmix.abstol'/'control$pnvmix.reltol' has not been reached.
##' @param ... additional arguments passed to the underlying mixing distribution  
##' @return numeric vector with the computed probabilities and attributes "error"
##'         (error estimate of the RQMC estimator) and "numiter"
##'         (number of iterations)
##' @author Erik Hintz 
corgnvmix <- function(scale, qmix, method = c("kendall", "spearman"),
                      groupings = 1:2, ellip.kendall = FALSE, control = list(), 
                      verbose = TRUE, ...)
{
   
   ## 1. Checks and set-up  ####################################################
   stopifnot(all(scale >= -1), all(scale <= 1), is.logical(ellip.kendall),
             is.logical(verbose))
   method <- match.arg(method) # note: 'method' used in different context below
   do.Kendall <- (method == "kendall")
   l.scale <- length( scale <- as.vector(scale) )
   ## Use approximate formula ktau = 2/pi * arcsin(scale) if required
   if(method == "kendall" & ellip.kendall){
      res <- 2/pi * asin(scale)
      #attr(res, "abs. error") <- rep(0, l.scale)
      #attr(res, "rel. error") <- rep(0, l.scale)
      #attr(res, "numiter") <- 0
      return(res)
   }
   ## Prepare variables for RQMC below 
   mix_list  <- get_mix_(qmix = qmix, groupings = groupings, 
                         callingfun = "dependencemeasures", ... )
   mix_      <- mix_list$mix_ # function(u) 
   numgroups <- length(unique(groupings)) # number of groups 
   control   <- get_set_param(control)
   ## Grab the following variables for readability
   method    <- control$method
   increment <- control$increment
   B         <- control$B
   ## Generate 'B' seeds for sobol(.., seed = seeds[b])
   seeds <- if(method == "sobol") sample(1:1e3, B) else NULL 
   ## Absolute or relative precision required?
   tol <- if(is.na(control$dependencemeasures.abstol)){
      do.reltol <- TRUE # use relative tolerance
      stopifnot(control$dependencemeasures.reltol > 0)
      control$dependencemeasures.reltol
   } else {
      do.reltol <- FALSE
      control$dependencemeasures.abstol
   }
   rqmc.estimates   <- matrix(0, ncol = l.scale, nrow = B)
   CI.factor.sqrt.B <- control$CI.factor / sqrt(B)
   ## Additional variables needed if the increment chosen is "doubling"
   if(increment == "doubling") {
      if(method == "sobol") useskip <- 0
      denom <- 1
   }
   ## Initialize max.error to > tol to enter while loop
   max.error <- tol + 42
   total.fun.evals <- 0
   numiter <- 0
   current.n <- control$fun.eval[1]
   dim <- if(do.Kendall) 2 else 3 # dimension of the point-sets used
   const <- if(do.Kendall) 2/pi else 6/pi 
   
   ## 2. Actual computation ####################################################
   ## while() runs until precision 'tol' is reached or the number of function
   ## evaluations exceed fun.eval[2] or the number of iterations exceed
   ## control$max.iter.rqmc. In each iteration, B RQMC estimates of
   ## of the dependence measure are computed; if 'scale' is a vector
   ## the same mixing realizations are used for all elements in 'scale'
   
   while(max.error > tol & total.fun.evals < control$fun.eval[2] &
         numiter < control$max.iter.rqmc)
   {
      
      ## Get B RQCM estimates
      for(b in 1:B)
      {
         ## 2.1 Get the point set ###########################################
         U <- switch(method,
                     "sobol" = {
                        if(increment == "doubling") {
                           qrng::sobol(n = current.n, d = dim,
                                       randomize = "digital.shift",
                                       seed = seeds[b],
                                       skip = (useskip * current.n))
                        } else {
                           qrng::sobol(n = current.n, d = dim,
                                       randomize = "digital.shift",
                                       seed = seeds[b],
                                       skip = (numiter * current.n))
                        }
                     },
                     "ghalton" = {
                        qrng::ghalton(n = current.n, d = dim,
                                      method = "generalized")
                     },
                     "PRNG" = {
                        matrix(runif(dim * current.n ), ncol = dim)
                     })
         
         ## 2.2 Evaluate the integrand at the (next) point set #################
         
         ## Obtain realizations of mix_ 
         mix.sample <- lapply(1:dim, function(i) {
            sample <- mix_(U[, i])
            if(!is.matrix(sample)) sample <- cbind(sample) # 1 or 2 col matrix
            sample[, groupings] # (current.n, 2) matrix 
         })
         ## Compute the terms involving 'mix.sample' (without 'scale') in asin()
         asinarg <- if(do.Kendall){
            (apply(mix.sample[[1]], 1, prod) + apply(mix.sample[[2]], 1, prod)) /
               sqrt( (mix.sample[[1]][, 1]^2 + mix.sample[[2]][, 1]^2) *
                        (mix.sample[[1]][, 2]^2 + mix.sample[[2]][, 2]^2))
         } else {
            sqrt(apply(mix.sample[[1]], 1, prod) / 
                    ( (mix.sample[[1]][, 1] + mix.sample[[2]][, 1]) * 
                         (mix.sample[[1]][, 2] + mix.sample[[3]][, 2]) ))
         } # current.n - vector 
         ## Compute next estimate 
         next.estimate <- const * .colMeans(asin(tcrossprod(asinarg, scale)), 
                                            m = current.n,  n = l.scale)
         ## 2.3 Update RQMC estimates ##########################################
         rqmc.estimates[b, ] <-
            if(increment == "doubling") {
               ## In this case both, rqmc.estimates[b] and
               ## next.estimate depend on n.current points
               (rqmc.estimates[b, ] + next.estimate) / denom
            } else {
               ## In this case, rqmc.estimates[b] depends on
               ## numiter * n.current points whereas next.estimate
               ## depends on n.current points
               (numiter * rqmc.estimates[b, ] + next.estimate) / (numiter + 1)
            }
      } # end for(b in 1:B)
      
      ## Update of various variables
      ## Number of function evaluations
      total.fun.evals <- total.fun.evals + B * current.n
      if(increment == "doubling") {
         ## Change 'denom' and 'useksip' (exactly once, in the first iteration)
         if(numiter == 0) {
            denom <- 2
            useskip <- 1
         } else {
            ## Increase sample size. This is done in all iterations
            ## except for the first two
            current.n <- 2 * current.n
         }
      }
      ## Update error depending on 'do.reltol'
      error <- if(!do.reltol) { # absolute error
         CI.factor.sqrt.B * apply(rqmc.estimates, 2, sd)
      } else { # relative error
         CI.factor.sqrt.B * apply(rqmc.estimates, 2, sd)/
            .colMeans(rqmc.estimates, m = B, n = l.scale)
      }
      max.error <- max(error)
      numiter <- numiter + 1 # update counter
   } # while ()
   
   ## 3. Finalize and return ###################################################
   res <- .colMeans(rqmc.estimates, m = B, n = l.scale)
   ## Handle warnings
   reached <- (error <= tol)
   if(any(!reached) & verbose > 0) {
      ii <- which(!reached)
      if(verbose == 1) {
         strng <- if(length(ii) > 6) {
            paste0(paste(head(ii), collapse = ", "), ",...")
         } else {
            paste(ii, collapse = ", ")
         }
         warning("Tolerance not reached for entries ",strng," of 'scale'; consider increasing 'fun.eval[2]' and 'max.iter.rqmc' in the 'control' argument.")
      } else {
         for(i in 1:length(ii)) {
            warning(sprintf("Tolerance not reached for entries %d of 'scale'; consider increasing 'fun.eval[2]' and 'max.iter.rqmc' in the 'control' argument", ii[i]))
         }
      }
   }
   ## Compute absolute and relative error for attributes
   abserror <- if(do.reltol){
      relerror <- error
      error * res
   } else {
      relerror <- error / res # 'error' is absolute error
      error
   }
   ## Return
   attr(res, "abs. error") <- abserror
   attr(res, "rel. error") <- relerror
   attr(res, "numiter") <- numiter
   res
}



##' @title Tail dependence coefficient for the grouped/ungrouped t copula 
##' 
##' @param df scalar or 2-vector giving degrees-of-freedoms of the copula 
##' @param scale n-vector giving the 'scale' parameters
##' @param control list(), see get_set_param()
##' @param verbose logical (or integer: 0 = FALSE, 1 = TRUE, 2 = more output)
##'        indicating whether a warning is given if the required precision
##'        'control$pnvmix.abstol'/'control$pnvmix.reltol' has not been reached.
##' @return numeric vector with the computed probabilities and attributes "error"
##'         (error estimate of the RQMC estimator) and "numiter"
##'         (number of iterations)
##' @author Erik Hintz 

lambda_gStudent <- function(df, scale, control = list(), verbose = TRUE)
{
   
   ## 1. Checks and set-up  ####################################################
   stopifnot(all(df > 0), all(scale > -1), all(scale < 1), is.logical(verbose))
   ## 'df' can have length 1 (ungrouped) or length 2 (grouped)
   l.df <- length( df <- as.vector(df) )
   stopifnot(1 <= l.df, l.df <= 2)
   l.scale <- length(scale <- as.vector(scale)) # 'scale' can be a vector
   
   ## In ungrouped case, closed formula for lambda is available
   if(l.df == 1){
      res <- 2 * pt( -sqrt( (df + 1)*(1 - scale) / (1 + scale)), df = df + 1)
      attr(res, "abs. error") <- rep(0, l.scale)
      attr(res, "rel. error") <- rep(0, l.scale)
      return(res)
   }
   control <- get_set_param(control)
   ## Absolute or relative tolerance?
   tol <- if(is.na(control$dependencemeasures.abstol)){
      do.reltol <- TRUE # use relative tolerance
      stopifnot(control$dependencemeasures.reltol > 0)
      control$dependencemeasures.reltol
   } else {
      do.reltol <- FALSE
      control$dependencemeasures.abstol
   }
   ## Grab method and prepare variables for RQMC
   method <- control$method
   B      <- control$B
   increment <- control$increment
   if(method == "sobol") seeds <- sample(1:(1e5*B), B) # B seeds for 'sobol()'
   rqmc.estimates   <- matrix(0, ncol = l.scale, nrow = B)
   CI.factor.sqrt.B <- control$CI.factor / sqrt(B)
   error <- rep(NA, l.scale) # store error for all 'scale'
   ## Additional variables needed if the increment chosen is "doubling"
   if(increment == "doubling") {
      if(method == "sobol") useskip <- 0
      denom <- 1
   }
   ## Initialize max.error to > tol to enter while loop
   max.error <- tol + 42
   total.fun.evals <- 0
   numiter <- 0
   current.n <- control$fun.eval[1]
   ## Constant needed repeatedly
   B_df <-
      c((2^(df[2]/2-df[1]/2) * gamma((1+df[2])/2) / gamma((1+df[1])/2))^(1/df[2]),
        (2^(df[1]/2-df[2]/2) * gamma((1+df[1])/2) / gamma((1+df[2])/2))^(1/df[1]))
   
   ## 2. Actual computation ####################################################
   ## while() runs until precision 'tol' is reached or the number of function
   ## evaluations exceed fun.eval[2] or the number of iterations exceed
   ## control$max.iter.rqmc. In each iteration, B RQMC estimates of
   ## of the tail dependence coefficient lambda are computed; if 'scale' is a vector
   ## the same mixing realizations are used for all elements in 'scale'
   
   while(max.error > tol & total.fun.evals < control$fun.eval[2] &
         numiter < control$max.iter.rqmc)
   {
      
      ## Get B RQCM estimates
      for(b in 1:B)
      {
         ## 2.1 Get the point set ###########################################
         U <- switch(method,
                     "sobol" = {
                        if(increment == "doubling") {
                           qrng::sobol(n = current.n, d = 1,
                                       randomize = "digital.shift",
                                       seed = seeds[b],
                                       skip = (useskip * current.n))
                        } else {
                           qrng::sobol(n = current.n, d = 1,
                                       randomize = "digital.shift",
                                       seed = seeds[b],
                                       skip = (numiter * current.n))
                        }
                     },
                     "ghalton" = {
                        qrng::ghalton(n = current.n, d = 1,
                                      method = "generalized")
                     },
                     "PRNG" = {
                        matrix(runif( current.n ), ncol = 1)
                     })
         ## 2.2 Evaluate the integrand at the (next) point set #################
         ## Realizations chi^2(df[1]+1) and chi^2(df[2]+1) via qgamma():
         chisq_sample <-
            cbind(qgamma(U, shape = (df[1] + 1)/2, scale = 2),
                  qgamma(U, shape = (df[2] + 1)/2, scale = 2))
         ## Compute next estimate based on 'chisq_sample'
         next.estimate <-
            colMeans(sapply(1:l.scale, function(i){
               pnorm(- (B_df[1] * chisq_sample[, 1]^(df[1]/(2*df[2])) -
                           scale[i]*sqrt(chisq_sample[, 1]))*(1-scale[i]^2)^(-1/2)) +
                  pnorm(- (B_df[2] * chisq_sample[, 2]^(df[2]/(2*df[1])) -
                              scale[i]*sqrt(chisq_sample[, 2])) *(1-scale[i]^2)^(-1/2))
            })) # l.scale - vector
         
         ## 2.3 Update RQMC estimates ##########################################
         rqmc.estimates[b, ] <-
            if(increment == "doubling") {
               ## In this case both, rqmc.estimates[b] and
               ## next.estimate depend on n.current points
               (rqmc.estimates[b, ] + next.estimate) / denom
            } else {
               ## In this case, rqmc.estimates[b] depends on
               ## numiter * n.current points whereas next.estimate
               ## depends on n.current points
               (numiter * rqmc.estimates[b, ] + next.estimate) / (numiter + 1)
            }
      } # end for(b in 1:B)
      
      ## Update of various variables
      ## Number of function evaluations
      total.fun.evals <- total.fun.evals + B * current.n
      if(increment == "doubling") {
         ## Change 'denom' and 'useksip' (exactly once, in the first iteration)
         if(numiter == 0) {
            denom <- 2
            useskip <- 1
         } else {
            ## Increase sample size. This is done in all iterations
            ## except for the first two
            current.n <- 2 * current.n
         }
      }
      ## Update error depending on 'do.reltol'
      error <- if(!do.reltol) { # absolute error
         CI.factor.sqrt.B * apply(rqmc.estimates, 2, sd)
      } else { # relative error
         CI.factor.sqrt.B * apply(rqmc.estimates, 2, sd)/
            .colMeans(rqmc.estimates, m = B, n = l.scale)
      }
      max.error <- max(error)
      numiter <- numiter + 1 # update counter
   } # while ()
   
   ## 3. Finalize and return ###################################################
   res <- .colMeans(rqmc.estimates, m = B, n = l.scale)
   ## Handle warnings
   reached <- (error <= tol)
   if(any(!reached) & verbose > 0) {
      ii <- which(!reached)
      if(verbose == 1) {
         strng <- if(length(ii) > 6) {
            paste0(paste(head(ii), collapse = ", "), ",...")
         } else {
            paste(ii, collapse = ", ")
         }
         warning("Tolerance not reached for entries ",strng," of 'scale'; consider increasing 'fun.eval[2]' and 'max.iter.rqmc' in the 'control' argument.")
      } else {
         for(i in 1:length(ii)) {
            warning(sprintf("Tolerance not reached for entries %d of 'scale'; consider increasing 'fun.eval[2]' and 'max.iter.rqmc' in the 'control' argument", ii[i]))
         }
      }
   }
   ## Compute absolute and relative error for attributes
   abserror <- if(do.reltol){
      relerror <- error
      error * res
   } else {
      relerror <- error / res # 'error' is absolute error
      error
   }
   ## Return
   attr(res, "abs. error") <- abserror
   attr(res, "rel. error") <- relerror
   attr(res, "numiter") <- numiter
   res
}
