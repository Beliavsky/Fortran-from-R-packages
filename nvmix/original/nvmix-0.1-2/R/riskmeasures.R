### Risk measures ###############################################################

##' @title Estimating value-at-risk for univariate normal variance mixtures
##' @param level vector of confidence levels
##' @param qmix see ?pnvmix()
##' @param loc see ?pnvmix()
##' @param scale see ?pnvmix()
##' @param control see ?get_set_param()
##' @param verbose logical if warnings should be thrown
##' @param ... see ?pnvmix()
##' @return vector of expected shortfall estimates with attributes 'error'
##'         and 'numiter'
##' @author Erik Hintz, Marius Hofert and Christiane Lemieux
ES_nvmix <- function(level, qmix, loc = 0, scale = 1, control = list(),
                    verbose = TRUE, ...){
   ## 1. Checks and variable declarations ######################################
   stopifnot(scale > 0)
   if(!is.vector(level)) level <- as.vector(level)
   if(any(level >= 1) | any(level <= 0)) stop("All levels in 'level' must be in (0,1)")
   ## Deal with algorithm parameters, see also ?get_set_param()
   control <- get_set_param(control)
   ## Define the quantile function of the mixing variable #######################
   mix_list      <- get_mix_(qmix = qmix, callingfun = "ESnvmix", ... )
   qW            <- mix_list[[1]] # function(u)
   special.mix   <- mix_list[[2]]
   ## In case of normal and t, the expected shortfall is known
   if(!is.na(special.mix) && !(special.mix == "pareto")){
      res <- switch(special.mix,
                    "inverse.gamma" = {
                       df <- mix_list$param
                       loc + sqrt(scale)*dt(qt(level, df = df), df = df)/(1-level)*
                          ((df + qt(level, df = df)^2)/(df-1))},
                    "constant" = {
                       loc + sqrt(scale)*dnorm(qnorm(level))/(1-level)
                    })
      numiter <- 0
      relerror <- rep(0, length(level))
      abserror <- rep(0, length(level))
   } else {
      ## Otherwise use RQMC to estimate the expected shortfall
      ## Estimate/compute VaR_alpha first
      VaRs <- qnvmix(level, qmix = qmix, control = control, ...)
      ## Initialize various quanitities
      total.fun.evals <- 0
      numiter <- 0
      method  <- control$method
      increment <- control$increment
      B <- control$B
      current.n <- control$fun.eval[1]
      ## Absolte/relative precision?
      tol <- if(is.na(control$riskmeasures.abstol)) {  # use relative error
         stopifnot(control$riskmeasures.reltol > 0)
         do.reltol <- TRUE
         control$riskmeasures.reltol
      } else {
         do.reltol <- FALSE # use absolute error
         control$riskmeasures.abstol
      }
      ## Additional variables needed if the increment chosen is "doubling"
      if(increment == "doubling") {
         if(method == "sobol") useskip <- 0
         denom <- 1
      }
      ## Sample 'B' seeds for 'sobol(..., seed = seeds_[b])' to get the same shifts
      if(method == "sobol") seeds_ <- sample(1:(1e5*B), B) # B seeds for 'sobol()'
      ## Matrix to store RQMC estimates for all levels in the vector 'level'
      rqmc.estimates <- matrix(0, ncol = length(level), nrow = B)
      ## Will be needed a lot:
      CI.factor.sqrt.B <- control$CI.factor / sqrt(B)
      sqrt.scale <- sqrt(scale)
      sqrt.two.pi <- sqrt(2*pi)
      ## Initialize error to > tol to enter while loop
      error <- tol + 42
      ## 2. Actual computation #################################################
      ## while() runs until precision 'tol' is reached or the number of function
      ## evaluations exceed fun.eval[2]. In each iteration, B RQMC estimates of
      ## the expected shortfall are computed; if 'level' is a vector,
      ## the same mixing realizations are used for all levels
      while(max(error) > tol & total.fun.evals < control$fun.eval[2] &
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
                                          seed = seeds_[b],
                                          skip = (useskip * current.n))
                           } else {
                              qrng::sobol(n = current.n, d = 1,
                                          randomize = "digital.shift",
                                          seed = seeds_[b],
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
            ## Realizations of sqrt(W)
            sqrt.mixings <- sqrt(qW(U))

            ## 2.2 Evaluate the integrand at the (next) point set ##############
            next.estimate <- sqrt.scale *
               .colMeans(sqrt.mixings*exp(-tcrossprod(1/sqrt.mixings, VaRs)^2/2)/
                            sqrt.two.pi,
                         m = current.n, n = length(level))/(1-level)

            ## 2.3 Update RQMC estimates #######################################
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
               .colMeans(rqmc.estimates, m = B, n = length(level))
         }
         numiter <- numiter + 1 # update counter
      } # while ()
      res <- loc + .colMeans(rqmc.estimates, m = B, n = length(level))
      ## Handle warnings
      reached <- (error <= tol)
      if(any(!reached) && verbose > 0) {
         ii <- which(!reached)
         if(verbose == 1) {
            strng <- if(length(ii) > 6) {
               paste0(paste(head(ii), collapse = ", "), ",...")
            } else {
               paste(ii, collapse = ", ")
            }
            warning("Tolerance not reached for entries ",strng," of confidence levels; consider increasing 'fun.eval[2]' and 'max.iter.rqmc' in the 'control' argument.")
         } else {
            for(i in 1:length(ii)) {
               warning(sprintf("Tolerance not reached for entries %d of confidence levels; consider increasing 'fun.eval[2]' and 'max.iter.rqmc' in the 'control' argument", ii[i]))
            }
         }
      }
      abserror <- if(do.reltol){
         relerror <- error
         error * res
      } else {
         relerror <- error / res # 'error' is absolute error
         error
      }
   } # else

   ## 3. Return ################################################################

   attr(res, "abs. error") <- abserror
   attr(res, "rel. error") <- relerror
   attr(res, "numiter") <- numiter
   res
}


##' @title Estimating value-at-risk for univariate normal variance mixtures
##' @param level vector of confidence levels
##' @param qmix see ?pnvmix()
##' @param loc see ?pnvmix()
##' @param scale see ?pnvmix()
##' @param control see ?get_set_param()
##' @param verbose logical if warnings should be thrown
##' @param ... see ?pnvmix()
##' @return vector of value at risk estimates
##' @author Erik Hintz, Marius Hofert and Christiane Lemieux
VaR_nvmix <- function(level, qmix, loc = 0, scale = 1, control = list(),
                     verbose = TRUE, ...){
   ## This is called by qnvmix(level, ...)
   loc + sqrt(scale) * quantile_(level, qmix = qmix, which = "nvmix1", d = 1,
                               control = control, verbose = verbose,
                               q.only = TRUE, stored.values = NULL, ...)
}



