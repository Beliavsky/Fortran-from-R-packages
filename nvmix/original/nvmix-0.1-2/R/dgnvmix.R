### dgnvmix() ###################################################################

##' @title Density of a Grouped  Normal Variance Mixture for restricted W
##' @param qW function of one variable specifying the quantile functions of W.
##' @param x n*d vector of evaluation points for the density 
##' @param factor.inv lower triangular part of the inverse of 'factor' 
##' @param lrdet log(sqrt(det(scale)))
##' @param u.left numeric in (0,1)
##' @param u.right numeric in (0,1), > u.left. Density will be estimated
##'         conditional on W being between its 'u.left' and 'u.right' quantile.
##' @param groupings see ?dgnvmix()         
##' @param max.iter.rqmc maximum number of iterations
##' @param return.all logical; if true, matrix (U, qW(U)) also returned.
##' @param control see ?get_set_param() 
##' @return List of three:
##'         $ldensities n-vector with computed log-density values
##'         $numiter numeric, number of iterations needed
##'         $error n-vector of error estimates for log-densities; either
##'         relative error or absolte error depending on is.na(control$dnvmix.reltol)
##'         $U_lh (N, n+1) matrix with columns 'u' and 'logh(u)' (for all rows in 'x')
##'         (only if return.all = TRUE)
##' @author Erik Hintz and Marius Hofert
densgmix_rqmc <- function(qW, x, factor.inv, lrdet = 0, u.left = 0, u.right = 1, 
                          groupings = 1:d, max.iter.rqmc = 10, return.all = FALSE,
                          control = list())
{
   ## 1. Checks and set-up  ####################################################
   stopifnot(is.function(qW)) # sanity check
   if (!is.matrix(x)) x <- rbind(x)
   d <- ncol(x)
   n <- nrow(x)
   numgroups <- length(unique(groupings)) # number of groups 
   ## Is there an integration region? 
   if(u.left == u.right){
      return(list(ldens = -Inf, error = 0, numiter = 0))
   }
   ## Define various quantities 
   control <- get_set_param(control)
   dblng   <- (control$increment == "doubling")
   B       <- control$B # number of randomizations
   current.n <- control$fun.eval[1] # initial sample size
   numiter   <- rep(0, n) # counter for the number of iterations
   max.numiter <- 0 # stores max(numiter) 
   total.fun.evals <- 0 # counter for the number of fun evaluations 
   ## Absolte/relative precision?
   if(is.na(control$dnvmix.reltol)) {
      ## Use absolute error
      tol <- control$dnvmix.abstol
      do.reltol <- FALSE
   } else {
      ## Use relative error (=> default with tol = 0.01) 
      tol <- control$dnvmix.reltol
      do.reltol <- TRUE
   }
   ## Store and create seeds if 'sobol' is used to get the same shifts later
   if(control$method == "sobol") {
      seeds_ <- sample(1:(1e3*B), B) # B seeds for 'sobol()'
   }
   ## Additional variables needed if the increment chosen is "dblng"
   if(dblng) {
      if(control$method == "sobol") useskip <- 0
      denom <- 1
   }
   ## Matrix to store RQMC estimates
   rqmc.estimates <- matrix(-Inf, ncol = n, nrow = B)
   ## Define trafo-function that maps u to (u.left, u.right) 
   trafo <- function(u) u.left + (u.right - u.left)*u
   ## Initialize 'max.error' to > tol so that we can enter the while loop
   max.error <- tol + 42
   notRchd <- 1:n # row indices in 'x' for which tol not reached 
   l.notRchd <- n 
   ## Matrix to store (u, W_1(u),...,W_k(u)) 
   if(return.all) {
      max.nrow <- if(dblng) current.n*B*2^(max.iter.rqmc-1) else 
         max.iter.rqmc*B*current.n
      U_lh <- matrix(NA, ncol = n+1, nrow = max.nrow)
      curr.lastrow <- 0 # counter row-index additional points are being inserted after
   }
   ## Some more constants needed multiple times
   lconst <- -d/2 * log(2*pi) - lrdet
   CI.factor.sqrt.B <- control$CI.factor / sqrt(B)
   ZERO <- .Machine$double.neg.eps # avoid evaluation at 0 < ZERO 
   firstiter <- TRUE # logical if we're in the first iteration 
   ZERO <- .Machine$double.neg.eps
   ## 2. Main loop #############################################################
   
   ## while() runs until precision abstol is reached or the number of function
   ## evaluations exceed fun.eval[2]. In each iteration, B RQMC estimates of
   ## the desired log-densities are calculated.
   while(l.notRchd > 0 & max.numiter < max.iter.rqmc &
         total.fun.evals < control$fun.eval[2])
   {
      ## In each randomization ...
      for(b in 1:B) {
         ## Get the point set (*not* sorted here!)
         U <- switch(
            control$method,
            "sobol" = {
               if(dblng) {
                  qrng::sobol(n = current.n, d = 1, randomize = "digital.shift", 
                              seed = seeds_[b], skip = (useskip * current.n))
               } else {
                  qrng::sobol(n = current.n, d = 1, randomize = "digital.shift", 
                              seed = seeds_[b], skip = (max.numiter * current.n))
               }
            },
            "ghalton" = {
               qrng::ghalton(n = current.n, d = 1, method = "generalized")
            },
            "PRNG" = {
               runif(current.n)
            }) 
         ## Obtain realizations of W_j, j = 1,..,numgroups 
         mixings <- pmax(qW((U <- trafo(U))), ZERO)
         ## Compute next estimate (and 'lhvals', if needed) in C 
         next.estimate <- if(return.all){
            lhvals <- .Call("eval_gdensmix_integrand_returnall",
                            x = as.double(as.vector(x[notRchd, ])),
                            mix = as.double(as.vector(mixings)),
                            groupings = as.integer(groupings),
                            factorinv = as.double(factor.inv),
                            d = as.integer(d),
                            N = as.integer(l.notRchd),
                            n = as.integer(current.n),
                            lconst = as.double(lconst))
            ## Compute next estimates via LogSumExp trick
            -log(current.n) + logsumexp(lhvals) # l.notRchd -vector
         } else { 
            -log(current.n) + .Call("eval_gdensmix_integrand_LSE",
                                    x = as.double(as.double(as.vector(x[notRchd, ]))),
                                    mix = as.double(as.vector(mixings)),
                                    groupings = as.integer(groupings),
                                    factorinv = as.double(factor.inv),
                                    d = as.integer(d),
                                    N = as.integer(l.notRchd),
                                    n = as.integer(current.n),
                                    lconst = as.double(lconst))
         }
            
         # 
         # if(!is.matrix(mixings)) mixings <- cbind(mixings) # can happen in ungrouped case
         # ## Transform to  1/sqrt(W_j), j=1,..,*d* 
         # mixings <- mixings[, groupings, drop = FALSE] # (current.n, d)
         # rt.mix.i <- t(1/sqrt(mixings)) # (d, current.n)  with1/sqrt(W_j), j=1,..,d  
         # ## Compute mahalanobis distances for each row in 'x' 
         # ## Done very often => use lapply and unlist rather than sapply(, simplify = TRUE)
         # mahasq <- array(unlist(lapply(notRchd, function(i){
         #    Dix <- rt.mix.i * matrix(x[i, ], ncol = current.n, nrow = d, byrow = FALSE)
         #    .colSums(Dix * (scale.inv %*% Dix), m = d, n = current.n)
         # })), dim = c(current.n, l.notRchd ))
         # # mahasq <- sapply(notRchd, function(i){
         # #    Dix <- rt.mix.i * matrix(x[i, ], ncol = current.n, nrow = d, byrow = FALSE)
         # #    .colSums(Dix * (scale.inv %*% Dix), m = d, n = current.n)
         # # }) # (current.n, l.notRchd ) matrix 
         # ## Compute values of log h(u) (i.e., the logarithmic integrand)
         # lhvals <- lconst - .rowSums(log(mixings), m = current.n, n = d)/2 - 
         #    mahasq/2 # (current.n, l.notRchd ) matrix 
         # ## Store if needed
         if(return.all) {
            U_lh[(curr.lastrow+1):(curr.lastrow+current.n), 1] <- U
            U_lh[(curr.lastrow+1):(curr.lastrow+current.n), 1 + notRchd] <- lhvals
            curr.lastrow <- curr.lastrow + current.n
         }

         ## Update RQMC estimates
         rqmc.estimates[b, notRchd] <-
            if(dblng) {
               ## In this case both, rqmc.estimates[b,] and
               ## next.estimate depend on n.current points
               .Call("logsumexp2",
                     a = as.double(rqmc.estimates[b, notRchd]),
                     b = as.double(next.estimate),
                     n = as.integer(l.notRchd )) - log(denom)
            } else {
               ## In this case, rqmc.estimates[b,] depends on
               ## numiter * current.n points whereas next.estimate
               ## depends on current.n points
               .Call("logsumexp2",
                     a = as.double(rqmc.estimates[b, notRchd] + log(max(numiter))),
                     b = as.double(next.estimate),
                     n = as.integer(l.notRchd )) - log(max.numiter + 1)
            }
      } # end for(b in 1:B)
      ## Update of various variables
      ## Double sample size and adjust denominator in averaging as well as useskip
      if(dblng) {
         ## Change denom and useksip (exactly once, in the first iteration)
         if(firstiter) {
            denom <- 2
            useskip <- 1
            firstiter <- FALSE 
         } else {
            ## Increase sample size n. This is done in all iterations
            ## except for the first two
            current.n <- 2 * current.n
         }
      }
      ## Total number of function evaluations
      total.fun.evals <- total.fun.evals + B * current.n
      numiter[notRchd] <- numiter[notRchd] + 1
      max.numiter <- max.numiter + 1
      ## Update error. The following is slightly faster than 'apply(..., 2, var)'
      ldens <- logsumexp(rqmc.estimates) - log(B) # performs better than .colMeans
      vars <- .colMeans((rqmc.estimates - rep(ldens, each = B))^2, B, n, 0)
      error <- if(!do.reltol) { # absolute error
         sqrt(vars)*CI.factor.sqrt.B
      } else { # relative error
         sqrt(vars)/abs(ldens)*CI.factor.sqrt.B
      }
      l.notRchd <- length(notRchd <- which(error > tol))
   }
   ## 3. Return ################################################################
   
   if(return.all) {
      list(ldensities = ldens, numiter = numiter, error = error,
           U_lh = U_lh[1:curr.lastrow,])
   } else {
      list(ldensities = ldens, numiter = numiter, error = error)
   }
}

##' @title Density of a Grouped Normal Variance Mixture
##' @param x (n, d)-matrix of evaluation points
##' @param qmix see ?pgnvmix()
##' @param loc see ?pgnvmix()
##' @param scale see ?pgnvmix()
##' @param factor if 'scale' not provided, scale = tcrossprod(factor)
##' @param factor.inv inverse of 'factor'
##' @param control list; see ?get_set_param()
##' @param log logical indicating whether the logarithmic density is to be computed
##' @param verbose logical indicating whether warnings shall be thrown.
##' @param ... additional arguments passed to the underlying mixing distribution
##' @return n-vector with computed density values and attributes 'error'
##'         (error estimate) and 'numiter' (number of while-loop iterations)
##' @note Computation of maha distances different in grouped case (=> NOT 
##'        ellitpical) => need 'scale.inv'         
##' @author Erik Hintz and Marius Hofert
dgnvmix <- function(x, groupings = 1:d, qmix, loc = rep(0, d), scale = diag(d),
                    factor = NULL, factor.inv = NULL, 
                    control = list(), log = FALSE, verbose = TRUE, ...)
{
   
   ## 1 Setup ##################################################################
   
   if(!is.matrix(x)) x <- rbind(x)
   d <- ncol(x) # dimension
   if(!is.matrix(scale)) scale <- as.matrix(scale)
   stopifnot(length(loc) == d, dim(scale) == c(d, d))
   verbose <- as.logical(verbose)
   numgroups <- length(unique(groupings))
   ## Deal with algorithm parameters
   control <- get_set_param(control)
   ## Deal with 'qmix' 
   mix_list      <- get_mix_(qmix = qmix, groupings = groupings, 
                             callingfun = "dgnvmix", ... ) 
   qW            <- mix_list[[1]] # function(u)
   special.mix   <- mix_list[[2]] # string or NA
   need.scale <- (numgroups == 1 & !is.na(special.mix)) 
   ## Build result object (log-density)
   lres <- rep(-Inf, (n <- nrow(x))) # n-vector of results
   abserror <- rep(NA, n)
   relerror <- rep(NA, n)
   notNA <- rowSums(is.na(x)) == 0
   lres[!notNA] <- NA
   x <- x[notNA,, drop = FALSE] # non-missing data (rows)
   n <- nrow(x) # update 
   numiter <- 0 # initialize counter 
   ## Deal with 'factor' and compute lrdet = log(sqrt(det(scale))) = log(det(factor)) 
   if(is.null(factor.inv)){
      if(is.null(factor)){
         factor <- t(chol(scale))
      } else if(!hasArg(scale) & need.scale){
         ## 'scale' not provided, but needed 
         scale <- tcrossprod(factor)
      }
      lrdet <- sum(log(diag(factor)))
      factor.inv <- solve(factor)
      factor.inv <- factor.inv[lower.tri(factor.inv, diag = TRUE)] 
   } else {
      stopifnot(all.equal(dim(factor.inv), c(d, d)))
      if(need.scale) scale <- tcrossprod(solve(factor.inv))
      lrdet <- -sum(log(diag(factor.inv))) # det(A) = 1/det(A^{-1}) 
      factor.inv <- factor.inv[lower.tri(factor.inv, diag = TRUE)] 
   }
   ## Remove 'loc' 
   if(any(loc != 0)) x <- sweep(x, 2L, loc)
   
   ## Absolte/relative precision?
   tol <- if(is.na(control$dnvmix.reltol)) { # if 'reltol=NA' use absolute precision
      do.reltol <- FALSE
      control$dnvmix.abstol
   } else { # otherwise use relative precision (default)
      do.reltol <- TRUE
      control$dnvmix.reltol
   }
   ## 2 Actual computation ####################################################
   
   numiter <- 0 # initialize counter 
   ## Deal with the different distributions
   if(numgroups == 1 & !is.na(special.mix)) { # case of a classical NVM dist'n 
      maha2 <- mahalanobis(x, center = FALSE, scale) # squared mahalanobis distances 
      lres[notNA] <- switch(
         special.mix,
         "inverse.gamma" = {
            df <- mix_list$param
            lgamma((df + d) / 2) - lgamma(df/2) - (d/2) * log(df*pi) - lrdet - 
               (df+d)/2 * log1p(maha2/df)
         },
         "constant" = {
            -(d/2) * log(2 * pi) - lrdet - maha2/2
         },
         "pareto" = {
            alpha <- mix_list$param
            log(alpha) - d/2*log(2*pi) - lrdet - (alpha+d/2)*log(maha2/2) +
               pgamma(maha2/2, scale = 1, shape = alpha+d/2, log.p = TRUE) + 
               lgamma(alpha+d/2)
         })
      if(!log) lres <- exp(lres) # already exponentiate
      relerror <- rep(0, length(maha2))
      abserror <- rep(0, length(maha2))
   } else {
      ## General case of a grouped or ungrouped normal variance mixture (=> RQMC)
      ## Apply non-adaptive RQMC on all inputs first
      rqmc.obj <- densgmix_rqmc(qW, x = x, factor.inv = factor.inv, u.left = 0,
                                u.right = 1, lrdet = lrdet, groupings = groupings, 
                                max.iter.rqmc = control$dnvmix.max.iter.rqmc.pilot,
                                return.all = TRUE, control = control)
      ## Extract results
      ldens   <- rqmc.obj$ldensities
      numiter <- rqmc.obj$numiter
      error   <- rqmc.obj$error
      if(any(error > tol)){
         ## Call adaptive procedure here
         ## Accuracy not reached for some inputs => use adaptive method there 
         if(control$dnvmix.doAdapt){
            lconst <- -d/2 * log(2*pi) - lrdet
            min.stratlength <- control$dnvmix.tol.stratlength
            ZERO <- .Machine$double.neg.eps # avoid evaluation at 0 < ZERO 
            ONE <- 1-.Machine$double.neg.eps # avoid evaluation at 1 > ONE  
            tol.bisec <- control$dnvmix.tol.bisec # vector of length 3
            notRchd <- which(error > tol)
            x. <- x[notRchd, , drop = FALSE]
            n.notRchd <- nrow(x.) 
            ## Grab realizations of (u, logh(u)) for all 'x[notRchd, ]'
            U_lh <- rqmc.obj$U_lh[, c(1, 1 + notRchd)] # first columns has 'u' 
            Us <- U_lh[, 1] # vector!
            ord <- order(Us)
            Us <- c(ZERO, Us[ord], ONE) # 'Us' are sorted and 'ZERO'/'ONE' included
            n.lh <- length(Us)
            ## Compute logh(u) for u = ZERO and u = ONE 
            mixings_ZO <- qW(c(ZERO, ONE))
            lhvals_ZO <- .Call("eval_gdensmix_integrand_returnall",
                               x = as.double(as.vector(x.)),
                               mix = as.double(as.vector(mixings_ZO)),
                               groupings = as.integer(groupings),
                               factorinv = as.double(factor.inv),
                               d = as.integer(d),
                               N = as.integer(n.notRchd),
                               n = as.integer(2),
                               lconst = as.double(lconst))
            # if(!is.matrix(mixings_ZO)) mixings_ZO <- cbind(mixings_ZO) 
            # mixings_ZO <- mixings_ZO[, groupings, drop = FALSE] # (2, d)
            # rt.mix.i_ZO <- t(1/sqrt(mixings_ZO)) # (d, 2) with 1/sqrt(W_j), j=1,..,d  
            # mahasq_ZO <- sapply(1:n.notRchd, function(i){ 
            #    Dix <- rt.mix.i_ZO * matrix(x.[i, ], ncol = 2, nrow = d, byrow = FALSE)
            #    .colSums(Dix * (scale.inv %*% Dix), m = d, n = 2)
            # }) # (2, n.notRchd) matrix: mahalanobis distances for each row in 'x.' 
            # lhvals_ZO <- -d/2 * log(2*pi) - lrdet - 
            #    .rowSums(log(mixings_ZO), m = 2, n = d)/2 - mahasq_ZO/2 # (2, n.notRchd) matrix 
            ## Store all log h(u) (same order as 'Us' above)
            lhvals <- rbind(lhvals_ZO[1, ], U_lh[ord, -1, drop = FALSE], 
                            lhvals_ZO[2, ]) # (n.lh, n.notRchd) matrix 
            stopifnot( nrow(lhvals) == n.lh ) # check
            ## Result objects
            ldens_adapt    <- rep(NA, n.notRchd)
            error_adapt    <- rep(NA, n.notRchd)
            numiters_adapt <- rep(NA, n.notRchd)
            ## Go through all rows in 'x.' (= all columns of 'lhvals') to find 
            ## limits for the adaptive procedure
            for(i in 1:n.notRchd){
               ## Initialize 
               u.left <- NA
               u.right <- NA
               l.max <- lhvals[ (ind.max <- which.max(lhvals[, i])), i] # *observed* maximum
               u.max <- Us[ind.max] 
               ## Maximum at left/right endpoint?
               if(ind.max == 1){
                  u.left <- 0 # maximum at left endpoint (or close to it)
               } else if(ind.max == n.lh){
                  u.right <- 1 # maximum at the right endpoint (or close to it)
               } 
               ## Tolerance above which RQMC is used 
               l.tol.int.lower <- l.max - control$dnvmix.order.lower * log(10)
               l.tol.int.lower <- l.max - 15 * log(10)
               ## Find candidates for 'u.left' and 'u.right' 
               candid.left <- rep(NA, 2)
               candid.right <- rep(NA, 2)
               if(any(lhvals[, i] > l.tol.int.lower)){
                  ## Indices of 'u's so corresponding log h(u) >  l.tol.int.lower
                  ind.gr <- which(lhvals[, i] > l.tol.int.lower) # length >= 1
                  num.ind.gr <- length(ind.gr)
                  ## Candidate interval for 'u.left' 
                  candid.left  <- c(max(0, Us[ind.gr[1] - 1]), # < thshold
                                    min(u.max, Us[ind.gr[1]])) # > thshold
                  ## Note: If ind.gr[1] = 1 => U[ind.gr[1] - 1] = numeric(0)
                  ## => max(0, U[ind.gr[1] - 1]) = 0
                  ## Candidate interval for 'u.right'
                  last.ind.gr <- tail(ind.gr, 1) 
                  candid.right <- if(last.ind.gr == n.lh){
                     ## Largest 'u' in 'Us' satisfies log h(u) > l.tol.int.lower
                     c(Us[n.lh], 1)
                  } else {
                     c(Us[last.ind.gr], Us[last.ind.gr+1])
                  }
               } else {
                  ## No obs > threshold 
                  if(!is.na(u.right)){
                     ## 'u.right' was set above => set 'u.left' if not set yet 
                     if(is.na(u.left)){
                        u.left <- if(u.right >= ONE) 1 - min.stratlength else u.max - 
                           min.stratlength
                     }
                  } else {
                     ## 'u.right' was not set
                     if(!is.na(u.left)){
                        ## But 'u.left' was 
                        u.right <- if(u.left == 0) min.stratlength else u.max + 
                           min.stratlength
                     } else {
                        ## Neither 'u.left' nor 'u.right' was set
                        u.left <- u.max - min.stratlength
                        u.right <- u.max + min.stratlength
                     } 
                  }
               }
               uLuR <- c(u.left, u.right) # potentially NA 
               ## Find 'u.left' and 'u.right' if not already set
               candids <- rbind(candid.left, candid.right)
               for(k in 1:2) {
                  curr.candid <- candids[k, ]
                  if(!is.na(uLuR[k])) next # 'u.left' or 'u.right' already set
                  uLuR[k] <- if(diff(curr.candid) <= tol.bisec[1]) {
                     curr.candid[2]
                  } else {
                     ## Use bisection 
                     convd <- FALSE
                     iter.bisec <- 0
                     while(!convd && iter.bisec < control$dnvmix.max.iter.bisec) {
                        iter.bisec <- iter.bisec + 1
                        ## Next point to check
                        u.next <- mean(curr.candid)
                        ## Compute log h(u.next)
                        mixings.next <- qW(u.next)
                        l.int.next <- .Call("eval_gdensmix_integrand_returnall",
                                           x = as.double(as.vector(x.[i, ])),
                                           mix = as.double(as.vector(mixings.next)),
                                           groupings = as.integer(groupings),
                                           factorinv = as.double(factor.inv),
                                           d = as.integer(d),
                                           N = 1L,
                                           n = 1L,
                                           lconst = as.double(lconst))
                        # if(!is.matrix(mixings.next)) 
                        #    mixings.next <- rbind(mixings.next)
                        # mixings.next <- mixings.next[, groupings, drop = FALSE]
                        # Dix <- t(1/sqrt(mixings.next) * x.[i, ])
                        # mahasq.next <- sum(Dix * (scale.inv %*% Dix))
                        # l.int.next <- -d/2 * log(2*pi) - lrdet - 
                        #    sum(log(mixings.next))/2 - mahasq.next/2
                        diff <- l.int.next - l.tol.int.lower
                        ## Update 'curr.candid' depending on sign of 'diff' and check convergence:
                        if(k == 1){
                           if(diff > 0) curr.candid[2] <- u.next else curr.candid[1] <- u.next
                        } else {
                           if(diff > 0) curr.candid[1] <- u.next else curr.candid[2] <- u.next
                        }
                        convd <- 
                           (abs(diff) < tol.bisec[3]) || (diff(curr.candid) < tol.bisec[1])
                     }
                     u.next
                  }
               }
               ## For readability
               u.left  <- if(uLuR[1] <= ZERO) 0 else uLuR[1]
               u.right <- if(uLuR[2] >= ONE)  1 else uLuR[2]
               ## Integrate the two regions outside (u.left, u.right) via trapezoidal rules
               ## ... (0, u.left):
               ldens.left <- if(u.left == 0) -Inf else {
                  ## 0 < u.left < 1 => Find obs in (0, u.left)
                  u_sml <- (Us <= u.left)
                  sum_u_sml <- sum(u_sml)
                  if(sum_u_sml > 1) {
                     ## Case 1: We have >1 observations in (0, u.left)
                     last_sml <- which(u_sml)[sum_u_sml]
                     weights <- abs(c(Us[1], diff(Us[1:last_sml])))
                     logsumexp(
                        as.matrix(log(weights) + (c(-Inf, lhvals[1:(last_sml-1), i]) + 
                                                     lhvals[1:last_sml, i])/2, ncol = 1))
                     
                  } else {
                     ## Case 2: No or only one observations in (0, u.left) 
                     log(u.left) + l.tol.int.lower - log(2)
                  }
               }
               ## ... (u.right, 1):
               ldens.right <- if(u.right == 1) -Inf else {
                  ## 0 < u.right < 1 => Find obs in (u.right, 1)
                  u_gtr <- (Us >= u.right)
                  sum_u_gtr <- sum(u_gtr)
                  if(sum_u_gtr > 1) {
                     ## Case 1: We have >1 observations in (u.right, 1)
                     first_gtr <- which(u_gtr)[1]
                     weights <- abs(c(Us[1], diff(Us[first_gtr:n.lh])))
                     logsumexp(
                        as.matrix(log(weights) + (c(lhvals[(first_gtr+1):n.lh, i], -Inf) + 
                                                     lhvals[first_gtr:n.lh, i])/2, ncol = 1))
                  } else {
                     ## Case 2: No or only one observations in (u.right, 1)
                     log1p(-u.right) + l.tol.int.lower - log(2)
                  }
               }
               ## Integrate from 'u.left' to 'u.right' via RQMC  
               ldens.stratum.obj <- 
                  densgmix_rqmc(qW, x = x.[i, , drop = FALSE], lrdet = lrdet,
                                factor.inv = factor.inv,
                                u.left = u.left, u.right = u.right, groupings = groupings,
                                max.iter.rqmc = control$max.iter.rqmc - 
                                   control$dnvmix.max.iter.rqmc.pilot,
                                control = control)
               numiters_adapt[i] <- ldens.stratum.obj$numiter
               error_adapt[i] <- ldens.stratum.obj$error
               ldens_adapt[i] <- logsumexp(rbind(
                  ldens.left, ldens.right, ldens.stratum.obj$ldensities + 
                     log(u.right - u.left), deparse.level = 0))
            }
            ## End of adaptive procedure. Store results
            ldens[notRchd] <- ldens_adapt
            error[notRchd] <- error_adapt
            numiter[notRchd] <- numiter[notRchd] + numiters_adapt
            ## Handle warnings
            if(as.logical(verbose)) { 
               if(any(is.na(error))) {
                  ## At least one error is NA
                  warning("Estimation unreliable, corresponding error estimate NA")
               }
               whichNA <- which(is.na(error))
               if(any(error[setdiff(1:length(error), whichNA)] > tol)) # 'setdiff' needed if 'whichNA' is empty
                  warning("Tolerance not reached for all inputs; consider increasing 'max.iter.rqmc' in the 'control' argument.")
            }
         } else if (as.logical(verbose)){
            ## Adaptive method *not* used, print warning 
            warning("Tolerance not reached for all inputs; consider using the adaptive method by setting 'control$dnvmix.doAdapt' to 'TRUE'.")
         }
      }
      ## Compute absolute and relative error on log mu_hat ('ldens' has length(notNA))
      abserror[notNA] <- if(do.reltol){
         relerror[notNA] <- error
         relerror[notNA] * abs(ldens) 
      } else { # error is absolute error
         relerror[notNA] <- error / abs(ldens) 
         error 
      }
      ## Correct results and error if 'log = FALSE'
      if(!log){
         ldens <- exp(ldens)
         ## CI for mu: exp(logmu_hat +/- abserr(logmu_hat))) = (lower, upper)
         ## => compute max. error on mu_hat as max( (upper - mu), (mu - lower) ) 
         relerror[notNA] <- max( (exp(abserror[notNA]) - 1), (1 - exp(-abserror[notNA])) )
         abserror[notNA] <- ldens * relerror[notNA] # ldens already exponentiated 
      }
      ## Grab results, correct 'error' and 'lres' if 'log = FALSE'
      lres[notNA] <- ldens
   }
   
   ## 3. Return ################################################################
   
   ## Note that 'lres' was exponentiated already if necessary.
   attr(lres, "abs. error") <- abserror
   attr(lres, "rel. error") <- relerror
   attr(lres, "numiter") <- numiter
   lres
}
