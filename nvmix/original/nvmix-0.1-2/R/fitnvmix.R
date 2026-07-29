### fitnvmix() #################################################################

### Estimate weights E(1/W | X) ################################################

##' @title Estimate Weights for fitnvmix()
##' @param maha2.2 squared maha distances divided by 2 (length n)
##' @param qW specification of the quantile function of W as function(u, nu)
##' @param nu parameter (vector) 'nu' of W
##' @param lrdet log(sqrt(det(scale)))
##' @param d dimension
##' @param special.mix either NA or string. Currently supported are 'inverse.gamma'
##'         and 'pareto' in which cases analytical weights are calculated
##' @param control see ?fitnvmix()
##' @param verbose see ?fitnvmix()
##' @return List of three:
##'         $weights n-vector with computed log-density values
##'         $numiter numeric, number of iterations needed
##'         $error n-vector of error estimates for log-densities; either relative
##'         error or absolte error depending on is.na(control$dnvmix.reltol)
##'         $UsWs (B, n) matrix (U, qW(U)) where U are uniforms
##'         (only if return.all = TRUE)
##' @author Erik Hintz, Marius Hofert, Christiane Lemieux
##' @note Weights corresponds to delta_ki in the paper
weights_ <- function(maha2.2, qW, nu, lrdet, d, special.mix, control, verbose)
{
   verbose <- as.logical(verbose) # only logical needed here
   weights.warn.count <- 0 # record number of warnings when estimating weights
   if(!is.na(special.mix)) { # weights are known analytically
      weights <- switch(special.mix,
                        "inverse.gamma" = {
                           (nu + d) / (nu + maha2.2*2)
                        },
                        "pareto" = {
                           pgamma(1, shape = nu+d/2+1, scale = 1/maha2.2)/
                              pgamma(1, shape = nu+d/2, scale = 1/maha2.2)*
                              (nu + d/2)/maha2.2
                        })
      numiter <- 0
      error <- rep(0, length(maha2.2))
   } else { # weights need to be estimated
      ## Absolte/relative precision?
      if(is.na(control$weights.reltol)) {
         tol <- control$weights.abstol
         do.reltol <- FALSE
      } else {
         ## Use relative error
         tol <- control$weights.reltol
         do.reltol <- TRUE
      }
      ## lconst for the integrand
      lconst <- rep(-lrdet - d/2*log(2*pi), length(maha2.2))
      ## Call RQMC procedure to estimate weights non-adaptively
      rqmc.obj <-
         weights_rqmc(maha2.2, qW = qW, nu = nu, lconst = lconst, d = d,
                      max.iter.rqmc = control$dnvmix.max.iter.rqmc.pilot,
                      control = control, return.all = TRUE)
      ## Extract results
      weights <- rqmc.obj$weights
      numiter <- rep(rqmc.obj$numiter, length(maha2.2))
      error   <- rqmc.obj$error
      if(any(error > tol)) {
         ## Accuracy not reached for at least one 'maha2.2' value
         ## => Use adaptive procedure
         notRchd <- which(error > tol)
         qW. <- function(u) qW(u, nu = nu)
         ldens.obj <- densmix_adaptrqmc(qW., maha2.2 = maha2.2[notRchd],
                                        lconst = lconst[notRchd], d = d,
                                        UsWs = rqmc.obj$UsWs, control = control)
         lcond.obj <- densmix_adaptrqmc(qW., maha2.2 = maha2.2[notRchd],
                                        lconst = lconst[notRchd], d = d,
                                        k = d + 2, UsWs = rqmc.obj$UsWs,
                                        control = control)
         weights[notRchd]   <- exp(lcond.obj$ldensities - ldens.obj$ldensities)
         ## Which weights cannot be reliably estimated?
         which.errorNA     <- which(is.na(lcond.obj$error) | is.na(ldens.obj$error))
         if(any(which.errorNA)) {
            weights.warn.count <- 1
            if(verbose)
               warning("Some weights cannot be reliably estimated, corresponding error estimate NA")
            ## Check if 'weights' decreasing in 'maha2.2'
            n <- length(weights)
            for(i in 1:n) {
               if(i <= 2) {
                  next
               } else if(weights[i] <= weights[i-1]) {
                  next
               } else {
                  ## In this case, weights[i] > weights[i-1]
                  ## Case 1: Is there any weight beyond i which is smaller than weights[i-1]?
                  smallerthani <- which(weights[i:n] <= weights[i-1]) + i - 1
                  if(length(smallerthani) > 0) {
                     ## If that's the case, interpolate all weights between i
                     ## and the first one smaller than weights[i-1]
                     firstsmaller <- smallerthani[1]
                     slope <- (weights[firstsmaller] - weights[i-1]) /
                        (maha2.2[firstsmaller] - maha2.2[i-1])
                     weights[i:(firstsmaller-1)] <- weights[i-1] + slope*
                        (maha2.2[i:(firstsmaller-1)] - maha2.2[i-1])
                  } else {
                     ## Behavior of the function suggests log-log extrapolation:
                     slope <- log(weights[i-1] / weights[i-2]) /
                        log(maha2.2[i-1] / maha2.2[i-2])
                     weights[i:n] <- weights[i-1] *
                        exp(slope*log(maha2.2[i:n] / maha2.2[i-1]))
                  }
               }
            }
         }
      }
   }
   list(weights = weights, numiter = numiter, error = error,
        weights.warn.count = weights.warn.count)
}

##' @title Estimate weights for fitnvmix() via non-adaptive RQMC
##' @param maha2.2 squared maha distances divided by 2 (length n)
##' @param qW specification of the quantile function of W as function(u, nu)
##' @param nu parameter (vector) 'nu' of W
##' @param lconst see ?densmix_()
##' @param d dimension
##' @param control see ?fitnvmix()
##' @param max.iter.rqmc maximum number of RQMC iterations
##' @param return.all logical; if true, matrix (U, qW(U)) also returned.
##' @return List of three:
##'         $weights n-vector with computed log-density values
##'         $numiter numeric, number of iterations needed
##'         $error n-vector of error estimates for log-densities; either relative
##'         error or absolte error depending on is.na(control$dnvmix.reltol)
##'         $UsWs (B, n) matrix (U, qW(U)) where U are uniforms
##'         (only if return.all = TRUE)
##' @author Erik Hintz
##' @note corresponds to delta_ki in the paper
weights_rqmc <- function(maha2.2, qW, nu, lconst, d, max.iter.rqmc, control,
                         return.all)
{
   ## Define various quantites:
   dblng           <- TRUE
   B               <- control$B # number of randomizations
   n               <- length(maha2.2) # sample size
   current.n       <- control$fun.eval[1] #initial sample size
   numiter         <- 0 # counter for the number of iterations
   ZERO            <- .Machine$double.neg.eps
   total.fun.evals <- 0
   ## Absolte/relative precision?
   if(is.na(control$weights.reltol)) {
      ## Use absolute error
      tol <- control$weights.abstol
      do.reltol <- FALSE
   } else {
      ## Use relative error
      tol <- control$weights.reltol
      do.reltol <- TRUE
   }
   ## Store seed if 'sobol' is used to get the same shifts later:
   if(control$method == "sobol") {
      useskip <- 0 # to get correct shifts if method = "sobol"
      seeds_ <- sample(1:(1e3*B), B) # B seeds for 'sobol()'
   }
   denom <- 1
   ## Matrix to store RQMC estimates (for the log-weights)
   rqmc.estimates.lweights <- matrix(-Inf, ncol = n, nrow = B)
   ## Will be needed a lot
   CI.factor.sqrt.B <- control$CI.factor / sqrt(B)
   ## Initialize 'max.error' to > tol so that we can enter the while loop:
   max.error <- tol + 42
   ## Matrix to store U, W values => nrows = maximal number of funevals
   if(return.all) {
      max.nrow <- current.n*B*2^(max.iter.rqmc-1)
      UsWs <- matrix(NA, ncol = 2, nrow = max.nrow)
      curr.lastrow <- 0 # will count row-index after which additional points are being inserted
   }

   ## Main loop
   ## while() runs until precision 'tol' is reached or until the number of function
   ## evaluations/iterations exceeds fun.eval[2]/max.iter.rqmc.
   ## In each iteration, B RQMC estimates of the desired log-weights are calculated.
   while(max.error > tol && numiter < max.iter.rqmc &&
         total.fun.evals < control$fun.eval[2])
   {
      ## For each randomization
      for(b in 1:B) {
         ## Get the point set
         U <- sort(switch(control$method,
                          "sobol" = {
                             if(dblng) {
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
                             runif(current.n)
                          })) # sorted for later!
         ## Evaluate the integrand at the (next) point set
         W <- qW(U, nu = nu) # realizations of the mixing variable; sorted!
         ## Need to replace values < ZERO by ZERO. W is *sorted*, so check using
         ## loop instd of 'pmax' (more efficient)
         for(ind in 1:current.n) if(W[ind] < ZERO) W[ind] <- ZERO else break
         ## Update 'UsWs'
         if(return.all) {
            UsWs[(curr.lastrow + 1) : (curr.lastrow + current.n), ] <- cbind(U, W)
            curr.lastrow <- curr.lastrow + current.n
         }
         next.est.condexp <- .Call("eval_densmix_integrand",
                                   W          = as.double(W),
                                   maha2_2    = as.double(maha2.2),
                                   current_n  = as.integer(current.n),
                                   n          = as.integer(n),
                                   d          = as.integer(d),
                                   k          = as.integer(d + 2), # k=d+2 here!
                                   lconst     = as.double(lconst))
         next.est.ldens <- .Call("eval_densmix_integrand",
                                 W          = as.double(W),
                                 maha2_2    = as.double(maha2.2),
                                 current_n  = as.integer(current.n),
                                 n          = as.integer(n),
                                 d          = as.integer(d),
                                 k          = as.integer(d), # k=d here!
                                 lconst     = as.double(lconst))

         ## Update RQMC estimates
         rqmc.estimates.lweights[b, ] <-
            .Call("logsumexp2",
                  a = as.double(rqmc.estimates.lweights[b, ]),
                  b = as.double(next.est.condexp - next.est.ldens),
                  n = as.integer(n)) - log(denom)

      } # end for(b in 1:B)
      ## Update of various variables
      ## Double sample size and adjust denominator in averaging as well as useskip
      if(numiter == 0) {
         ## Change denom and useksip (exactly once, in the first iteration)
         denom <- 2
         useskip <- 1
      } else {
         ## Increase sample size n. This is done in all iterations
         ## except for the first two
         current.n <- 2 * current.n
      }
      ## Total number of function evaluations:
      total.fun.evals <- total.fun.evals + B * current.n
      numiter <- numiter + 1
      ## Update error for 'weights' (NOT log values!)
      rqmc.estimates.weights <- exp(rqmc.estimates.lweights)
      weights  <- colMeans(rqmc.estimates.weights)
      vars     <- .colMeans((rqmc.estimates.weights -
                                rep(weights, each = B))^2, B, n, 0)
      errors   <- if(!do.reltol) {
         sqrt(vars)*CI.factor.sqrt.B
      } else {
         sqrt(vars)/abs(weights)*CI.factor.sqrt.B
      }
      max.error <- max(errors)
   } # while()

   ## Return
   ret.obj <- if(return.all) {
      list(weights = weights, numiter = numiter, error = errors,
           UsWs = UsWs)
   } else {
      list(weights = weights, numiter = numiter, error = errors)
   }
   ret.obj
}


### Estimate 'nu' given 'loc' and 'scale' ######################################

##' @title Optimizer for univariate concave functions
##' @param fn function of one argument to be maximized
##' @param lower lower bound on position of the max
##' @param upper upper bound on position of the max
##' @param par.init initial value for optimization (by default (upper+lower)/2)
##' @param eps.bisec required length of starting interval found by
##'         'get_next_candid()' that is passed to optim
##' @param max.iter.bisec maximum number of calls to 'get_next_candid()'
##' @return list with elements 'par' and 'value' giving the maximum location
##'          and value of fn
##' @note   'optim_local()' essentially calles 'get_next_candid()' until an
##'          interval of length <= eps.bisec is found. This interval is then
##'          passed to 'optimize()'
##' @author Erik Hintz
optim1d_ <- function(fn, lower = 0.75, upper = 10, par.init = NULL,
                        eps.bisec = 0.25, max.iter.bisec = 20){
   ## 0. Handle 'par.init' if provided
   if(!is.null(par.init)){
      ## Make sure 'par.init' does not fall on boundary
      if(lower == par.init | upper == par.init) par.init <- (upper+lower)/2
      stopifnot(lower < par.init, par.init < upper) # sanity check
      ## Check fn() val's at 'par.init +/- eps.bisec'
      lower_ <- max(par.init - eps.bisec, lower)
      upper_ <- min(par.init + eps.bisec, upper)
      fnvals_ <- c(fn(lower_), fn(par.init), fn(upper_))
      ## Determine interval with maximum
      if(fnvals_[1] <= fnvals_[2] & fnvals_[2] >= fnvals_[2]){
         ## Max within par.init +/- eps.bisec
         candids_ <- c(lower_, par.init, upper_)
      } else if(fnvals_[1] <= fnvals_[2] & fnvals_[2] <= fnvals_[3]) {
         ## Max between par.init and (origninal) upper
         candids_ <- c(par.init, upper_, upper)
         fnvals_  <- c(fnvals_[2:3], fn(upper))
      } else if(fnvals_[1] >= fnvals_[2] & fnvals_[2] >= fnvals_[3]){
         ## Max between (original) lower and par.init
         candids_ <- c(lower, lower_, par.init)
         fnvals_  <- c(fn(lower), fnvals_[1:2])
      } else {
         ## Remaining case: fnvals[1] >= fnvals[2], fnvals[2] <= fnvals[3]
         ## But then fn() first decreasing, then increasing => can't be
         stop("Non-concavity detected")
      }
   } else {
      ## No 'par.init' provided
      par.init <- (upper+lower)/2
      candids_ <- c(lower, par.init, upper)
      fnvals_  <- sapply(1:3, function(i) fn(candids_[i]))
   }

   ## 1. Bisections to find starting interval
   ## Initialize length of candid interval to > eps.bisec to enter while() below
   l.candids <- candids_[3] - candids_[1] + eps.bisec
   iter.bisec <- 1
   while(l.candids > eps.bisec & iter.bisec < max.iter.bisec){
      next.candids.obj <-
         get_next_candid(fn, candid = candids_, fnvals = fnvals_)
      candids_ <- next.candids.obj$candid
      fnvals_  <- next.candids.obj$fnvals
      l.candids <- candids_[3] - candids_[1]
      iter.bisec <- iter.bisec + 1
   }

   ## 2. Optimization via 'optimize()'
   opt.obj <- optimize(fn, interval = c(candids_[1], candids_[3]), maximum = TRUE)

   ## 3. Return
   list(par = opt.obj$maximum, value = opt.obj$objective)
}


##' @title Bisection to find an interval with maximum value of a concave function
##' @param fn function of one argument to be maximized
##' @param candid 3-vector: 1st/3rd element give lower/upper bound
##'        on position of the max. 2nd element is an initial value
##' @param fnvals 3-vector: fn(candid)
##' @param verbose numeric or logical. 0: No warnings; 1: Warnings;
##'        2: Warnings + short tracing; 3: Warnings + complete tracing.
##' @return list with elements 'candid' and 'fnvals' where
##'         the new 'candid' as half as wide as the input 'candid'
##' @author Erik Hintz
get_next_candid <- function(fn, candid = c(1, 11/2, 10), fnvals = NULL){
   ## Compute function values at 'candid' if not supplied
   if(is.null(fnvals)) fnvals <- sapply(1:3, function(i) fn(candid[i]))
   par.next.l <- mean(candid[1:2]) # point in the 'left half'
   fn.next.l <- fn(par.next.l)
   if(fn.next.l > fnvals[2]){
      ## Found new interval (candid[1], candid[2])
      candid_ <- c(candid[1], par.next.l, candid[2])
      fnvals_ <- c(fnvals[1], fn.next.l, fnvals[2])
   } else {
      par.next.r <- mean(candid[2:3]) # point in the 'right half'
      fn.next.r <- fn(par.next.r)
      if(fn.next.r < fnvals[2]){
         ## Found new interval (par.next.l, par.next.r)
         candid_ <- c(par.next.l, candid[2], par.next.r)
         fnvals_ <- c(fn.next.l, fnvals[2], fn.next.r)
      } else {
         ## Found new interval (candid[2], candid[3])
         candid_ <- c(candid[2], par.next.r, candid[3])
         fnvals_ <- c(fnvals[2], fn.next.r, fnvals[3])
      }
   }
   ## Return
   list(candid = candid_, fnvals = fnvals_)
}



##' @title Estimate 'nu' given 'loc' and 'scale' by maximizing the Log-Likelihood
##' @param tx t(x) where x is as in ?fitnmvix()
##' @param qW specification of the quantile function of W as function(u, nu)
##' @param init.nu initial estimate of 'nu'
##' @param loc current estimate of 'loc'
##' @param scale current estimate of 'scale'
##' @param factor cholesky factor of the 'scale'; if not provided, it's calculated
##' @param mix.param.bounds see ?fitnvmix()
##' @param special.mix string specifying if W has a special distribution for which
##'        weights are known; currently, only 'inverse.gamma' and 'pareto' supported
##' @param control see ?fitnvmix()
##' @param control.optim passed to optim; see ?optim()
##' @param verbose see ?fitnvmix
##' @return list of two: $nu.est (scalar or vector of length 'init.nu'; MLE estimate of nu)
##'                      $max.ll (negative log-likelihood at nu.est)
##'                      $ll.counts (total number of calls to likelihood)
##'                      $opt.obj (object returned by underlying 'optim' call)
##' @author Erik Hintz, Marius Hofert, Christiane Lemieux
estim_nu <- function(tx, qW, init.nu, loc, scale, factor = NA, mix.param.bounds,
                     special.mix = NA, control, control.optim, verbose)
{
   ## Obtain 'factor' if not provided
   if(is.na(factor)) factor <- t(chol(scale))
   ll.counts <- 0 # counts number of calls to 'neg.log.likelihood.nu()'
   dnvmix.warn.count <- 0 # counts number of dnvmix warnings
   if(is.character(special.mix)) { # analytical density available
      stopifnot(special.mix == "inverse.gamma" || special.mix == "pareto")
      ## In this case, dnvmix() uses a closed formula for the density
      loglik <- function(nu) {
         ll.counts <<- ll.counts + 1 # update 'll.counts' in the parent environment
         if(any(nu < mix.param.bounds[, 1]) | any(nu > mix.param.bounds[, 2]))
            return(-Inf)
         ll <- sum(dnvmix(t(tx), qmix = special.mix, loc = loc, factor = factor,
                          nu = nu, log = TRUE, verbose = verbose))
         if(verbose >= 3) cat(".") # print dot after each call to likelihood
         ll # return
      }
   } else {
      ## Get various quantitites passed to 'densmix_()'
      z        <- forwardsolve(factor, tx - loc, transpose = FALSE)
      maha2.2  <- sort(colSums(z^2)/2)
      lrdet    <- sum(log(diag(factor)))
      d        <- ncol(factor)
      ## Generate a seed (=> results 'more monotone')
      seed_ <- sample(1:1e7, 1)
      ## Set up -loglikelihood as a function of 'nu'
      lconst <- rep(-lrdet - d/2*log(2*pi), length(maha2.2))
      loglik <- function(nu) {
         ll.counts <<- ll.counts + 1 # update 'll.counts' in parent environment
         if(any(nu < mix.param.bounds[, 1]) | any(nu > mix.param.bounds[, 2]))
            return(-Inf)
         set.seed(seed_) # reset seed => monotonicity
         qmix. <- function(u) qW(u, nu = nu) # function of u only
         ## Call 'densmix_()' which by default returns the log-density
         ldens.obj <- densmix_(qW = qmix., maha2.2 = maha2.2, lconst = lconst,
                               d = d, control = control, verbose = FALSE)
         ## Increase warning counter if necessary
         if(any(is.na(ldens.obj$relerror)))
            dnvmix.warn.count <<- dnvmix.warn.count + 1
         if(verbose >= 3) cat(".") # print dot after each call to likelihood
         ## Return -log-density
         sum(ldens.obj$ldensities)
      }
   }
   sign <- 1 # to report the correct sign of likelihood
   ## Optimize 'loglik' over 'nu'
   opt.obj <- if(dim(mix.param.bounds)[1] == 1){
      ## One dimensional parameter => manual search for *initial* nu, then optimize()
      optim1d_(loglik, lower = mix.param.bounds[1, 1],
               upper = mix.param.bounds[1, 2], par.init = init.nu)
   } else {
      ## Dimension > 1 => use optim with Nelder
      ## (works for non-differentiable functions)
      sign <- -1
      optim(init.nu, fn = function(nu) -loglik(nu), control = control.optim)
   }
   ## Return
   list(nu.est = opt.obj$par, max.ll = sign*opt.obj$value, ll.counts = ll.counts,
        opt.obj = opt.obj, dnvmix.warn.count = dnvmix.warn.count)
}


### Main function fitnvmix() ###################################################

##' @title Fitting Multivariate Normal Variance Mixtures
##' @param x (n, d) data matrix
##' @param qmix character string ("constant", "inverse.gamma", "pareto") or
##'        function(u, nu) interpreted as quantile function of the mixing rv W
##' @param mix.param.bounds either a vector of length two (in which case 'nu'
##'        is a scalar) or a matrix with two columns, where element in row i in
##'        1st/2nd column corresponds to lower/upper limits for component i of
##'        'nu'. All elements need to be finite, numeric values.
##'        (eg if nu/W~chi^2_nu, then eg 'mix.param.bounds = c(1, 10)'
##'        Note: The smaller the range, the better.
##' @param nu.init initial estimate for 'nu'; either NA in which case it is
##'        estimated or a vector of length = length of parameter vector 'nu'.
##'        If provided and close to MLE, can speed up 'fitnvmix' significantly
##' @param loc if provided, taken as the 'true' location vector
##' @param scale if provided, taken as the 'true' scale matrix
##' @param init.size.subsample if 'is.na(nu.init)', size of subsample of 'x'
##'        used to obtain initial estimate of nu.
##' @param size.subsample numeric, <= nrow(x). Number of rows of 'x' used in ECME
##'        iteration to optimize the log-likelihood. Defaults to n
##'        (all datapoints are used)
##' @param control list of algorithm specific parameters, see ?get_set_param
##'        and ?fitnvmix
##' @param verbose numeric or logical. 0: No warnings; 1: Warnings;
##'        2: Warnings + short tracing; 3: Warnings + complete tracing.
##' @return S3 object of class 'fitnvmix'; see below in 'class_fitnvmix()'
##' @author Erik Hintz, Marius Hofert, Christiane Lemieux
fitnvmix <- function(x, qmix, mix.param.bounds, nu.init = NA,
                     loc = NULL, scale = NULL,
                     init.size.subsample = min(n, 100), size.subsample = n,
                     control = list(), verbose = TRUE)
{
   ## 0 Setup ##################################################################
   call <- match.call() # for return
   if(!is.matrix(x))
      x <- rbind(x)
   ## Initialize various quantities
   control <- get_set_param(control)
   ## Prepare mixing variable
   mix_list      <- get_mix_(qmix = qmix, callingfun = "fitnvmix")
   qW            <- mix_list[[1]] # function(u, nu)
   special.mix   <- mix_list[[2]] # string or NA
   ## Check 'verbose' argument:
   if(is.logical(verbose)) {
      verbose <- as.integer(verbose)
   } else if(!(verbose %in% c(0, 1, 2, 3))) {
      stop("'verbose' has to be either logical or an integer between 0 and 3")
   }
   ## Check inputs, get dimensions
   notNA <- rowSums(is.na(x)) == 0
   x     <- x[notNA,, drop = FALSE] # non-missing data (rows)
   tx    <- t(x)
   n     <- nrow(x)
   d     <- ncol(x)
   ## Estimate 'loc' and 'scale'?
   do.loc   <- TRUE
   do.scale <- TRUE
   if(!is.null(loc)){
      stopifnot(length(loc <- as.vector(loc)) == d)
      do.loc <- FALSE
   }
   if(!is.null(scale)){
      if(!is.matrix(scale)) scale <- as.matrix(scale)
      stopifnot(dim(scale) == c(d, d))
      do.scale <- FALSE
   }
   ## Case of MVN: MLEs are sample mean and sample cov matrix
   is.mvn <- (is.character(special.mix) & special.mix == "constant")
   if(is.mvn) {
      loc.est <- if(do.loc) colMeans(x) else loc
      scale.est <- if(do.scale) as.matrix(nearPD(cov(x))$mat) else scale
      max.ll <- sum(dNorm(x, loc = loc.est, scale = scale.est, log = TRUE))
      return(class_fitnvmix(
         loc = loc.est, scale = scale.est, max.ll = max.ll, 
         data = x, 
         is.mvn = TRUE, do.loc = do.loc, do.scale = do.scale,
         call = call))
   }
   is.mvt <- # needed below to return 'df' instead of 'nu'
      (is.character(special.mix) & special.mix == "inverse.gamma")

   ## Use only sub-sample to estimate 'nu'?
   if(size.subsample < n) {
      sampled.ind <- sample(n, size.subsample)
      x.sub       <- x[sampled.ind,, drop = FALSE]
      tx.sub      <- tx[, sampled.ind, drop = FALSE]
   } else {
      x.sub  <- x
      tx.sub <- tx
   }
   ## Check parameter bounds on 'nu' and get 'mix.param.length'
   mix.param.length <- if(is.vector(mix.param.bounds)) {
      stopifnot(length(mix.param.bounds) == 2)
      mix.param.bounds <- matrix(mix.param.bounds, nrow = 1)
      1
   } else if(is.matrix(mix.param.bounds)) {
      stopifnot(dim(mix.param.bounds)[2] == 2)
      dim(mix.param.bounds)[1]
   } else stop("'mix.param.bounds' has to be either a vector of length 2 or a matrix with 2 columns")

   ## Matrix storing all 'nu' estimates with corresponding likelihoods
   nu.ests.ll <-
      matrix(NA, ncol = mix.param.length + 1,
             nrow = (nrow <- if(control$laststep.do.nu) control$ECME.maxiter + 2
                     else control$ECME.maxiter + 1))
   rownames(nu.ests.ll) <- if(control$laststep.do.nu) {
      c(paste0("Initial",
        sapply(1:control$ECME.maxiter,
               function(i) paste0("ECME-iteration ", i)),
        "Laststep"))
   } else {
      c(paste0("Initial (n0 = ", init.size.subsample, ")"),
        sapply(1:control$ECME.maxiter,
               function(i) paste0("ECME-iteration ", i)))
   }
   colnames(nu.ests.ll) <- c(sapply(1:mix.param.length,
                                 function(i) paste0("nu[", i, "]")),
                          "log-likelihood")
   current.iter.total <- 1

   ## Counters for warnings
   dnvmix.warn.count  <- 0
   weights.warn.count <- 0

   ## 1 Initial estimates for nu, loc, scale ##################################

   ## Unbiased estimator for 'loc' based on full sample
   loc.est <- if(do.loc) colMeans(x) else loc
   if(!do.scale) scale.est <- scale
   ## Sample covariance matrix based on full sample
   SCov <- as.matrix(nearPD(cov(x))$mat)
   ## Check if 'init.nu' was provided. If so, estimate 'scale' as (1/E(W))*SCov
   if(!is.na(nu.init)) {
      stopifnot(length(nu.init) == mix.param.length)
      if(verbose >= 2) cat("Step 1: Initial estimate for 'nu': Was provided")
      nu.est <- nu.init
      if(do.scale) scale.est <- 1/mean(qW(runif(1e4), nu.est))* SCov
      max.ll <- sum(dStudent(x, df = nu.est, scale = scale.est, log = TRUE))
   } else if(!do.scale){
      ## 'scale' was provided => only estimate 'nu'
      if(verbose >= 2)
         cat(paste0("Step 1: Initial estimate for 'nu' by optimizing log-likelihood over subsample of size ", init.size.subsample, " "))
      subsample <- x[sample(n, init.size.subsample),, drop = FALSE]
      opt.obj <-
         estim_nu(t(subsample), qW = qW, init.nu = rowMeans(mix.param.bounds),
                  loc = loc.est, scale = scale.est, mix.param.bounds = mix.param.bounds,
                  special.mix = special.mix, control = control,
                  control.optim = control$control.optim, verbose = verbose)
      nu.est <- opt.obj$nu.est
      max.ll <- opt.obj$max.ll
      dnvmix.warn.count <- opt.obj$dnvmix.warn.count + dnvmix.warn.count
   } else {
      ## Neither 'scale' nor 'nu' provided
      if(verbose >= 2)
         cat(paste0("Step 1: Initial estimate for 'nu' and 'scale' by optimizing log-likelihood over subsample of size ", init.size.subsample, " "))
      ll.counts <- 0 # counts number of calls to loglik function
      ## Generate and store seed
      seed_ <- sample(1:1e4, 1)
      ## Set up log-lik as function of 'nu' and 'c' where 'c' is scaling the
      ## argument 'scale = c*SCov'
      subsample <- x[sample(n, init.size.subsample),, drop = FALSE]
      loglik <- function(par) {
         ll.counts <<- ll.counts + 1 # update 'll.counts' in the parent environment
         ## Grab 'nu' and 'c' from 'par' for readability
         nu <- par[1:mix.param.length] # last element is 'c'
         c  <- par[mix.param.length+1]
         ## Return -Inf if 'par' outside range
         if(any(nu < mix.param.bounds[, 1]) | any(nu > mix.param.bounds[, 2]) |
            c < 0.1) return(-Inf)
         ## Define 'qmix_' as string (=> density known) or function
         qmix_ <- if(is.character(special.mix)) special.mix else qW
         set.seed(seed_) # same seed for different calls
         ll.obj <- dnvmix(subsample, qmix = qmix_, loc = loc.est,
                          scale = c * SCov, nu = nu, log = TRUE, verbose = FALSE)
         if(any(is.na(attr(ll.obj, "rel. error"))))
            dnvmix.warn.count <<- dnvmix.warn.count + 1
         if(verbose >= 3) cat(".") # print dot after each call to loglik
         sum(ll.obj) # return
      }
      ## Initial value for 'par' (was not provided)
      init.par <- c(rowMeans(mix.param.bounds),
                    1/mean(qW(runif(1e4), nu = rowMeans(mix.param.bounds))))
      ## Optimize log-lik
      opt.obj <- optim(init.par, fn = function(par) -loglik(par),
                       control = control$control.optim)
      ## Grab estimates for 'nu' and 'scale'
      nu.est    <- opt.obj$par[1:mix.param.length]
      scale.est <- opt.obj$par[mix.param.length + 1] * SCov
      max.ll    <- -opt.obj$value
      if(verbose >= 3)
         cat(paste0(" DONE (", ll.counts, " calls to likelihood needed)", '\n'))
   }
   ## Store additional values
   ## Matrix storing all 'nu' estimates with corresponding likelihoods
   nu.ests.ll[current.iter.total, ] <- c(nu.est, max.ll)
   current.iter.total <- current.iter.total + 1
   iter.converged <- control$ECME.maxiter + 2


   ## 2 ECME iteration ########################################################

   if(control$ECMEstep) {
      if(verbose == 2) cat(paste0('\n')) # if 'verbose==3' linebreak already happened
      if(verbose >= 2) cat(paste0("Step 2: ECME iteration.", '\n'))
      ## Initialize various quantities
      iter.ECME <- 0
      converged <- FALSE
      ## Main loop:
      while(iter.ECME < control$ECME.maxiter & !converged) {
         ## Print progress
         if(verbose >= 2) cat(paste0("  Iteration ",iter.ECME + 1, '\n'))
         ## 2.1  Update 'loc.est' and 'scale.est' while 'nu.est' held fixed ####
         if(do.loc | do.scale){ # otherwise can skip this step
            converged.locscale   <- FALSE
            iter.locscaleupdate  <- 1
            if(verbose >= 3)
               cat(paste0("    Estimating weights and updating 'loc' and 'scale'"))
            ## Inner loop (iterating over 'loc' and 'scale' with 'nu.est' held fixed)
            while(!converged.locscale &
                  iter.locscaleupdate < control$max.iter.locscaleupdate)
            {
               ## Get new maha distances (with current 'loc.est' and 'scale.est')
               factor <- t(chol(scale.est))
               lrdet <- sum(log(diag(factor)))
               z <- forwardsolve(factor, tx - loc.est, transpose = FALSE) # use the full sample!
               maha2.2.new <- colSums(z^2)/2
               order.maha2.2.new <- order(maha2.2.new)
               maha2.2.new <- maha2.2.new[order.maha2.2.new] # sorted increasingly
               if(iter.locscaleupdate == 1) {
                  ## Only in the first iteration do we approximate *all* weights by RQMC.
                  weights.obj <-
                     weights_(maha2.2.new, qW = qW, nu = nu.est, lrdet = lrdet, d = d,
                              special.mix = special.mix, control = control,
                              verbose = FALSE)
                  weights <- weights.obj$weights
                  weights.warn.count <-
                     weights.warn.count + weights.obj$weights.warn.count
                  weights.new <- weights[order(order.maha2.2.new)] # reorder
                  maha2.2     <- maha2.2.new # need to store maha-distances for interpolation
                  length.maha <- n # store length of 'maha2.2' and 'weights'
                  if(verbose >= 3)
                     cat(".") # print dot after estimation of weights.
               } else {
                  ## Linearly interpolate 'weights' to get new weights
                  weights.new          <- rep(NA, n)
                  curr.index           <- 1 # index to look for close values in 'maha2.2'
                  notInterpol          <- rep(NA, n)
                  notInterpolcounter   <- 1
                  for(ind in 1:n) {
                     curr.maha2.2 <- maha2.2.new[ind]
                     if(curr.maha2.2 < maha2.2[1] || curr.maha2.2 > maha2.2[length.maha]) {
                        ## Current maha too small or too large to use extrapolation
                        notInterpol[notInterpolcounter] <- ind
                        notInterpolcounter <- notInterpolcounter + 1
                     } else {
                        ## Start looking for close maha values in 'maha2.2'
                        found <- FALSE
                        while(!found && curr.index < length.maha) {
                           ## Found m1, m2 such that m1 <= curr.maha2.2 <= m2?
                           if(maha2.2[curr.index] <= curr.maha2.2 &
                              curr.maha2.2 <= maha2.2[curr.index+1]) {
                              found <- TRUE
                              ## Now check if we can interpolate (ie rel.error small)
                              if(abs(weights[curr.index+1] - weights[curr.index])/
                                 weights[curr.index+1] < control$weights.interpol.reltol)
                              {
                                 weights.new[ind] <- weights[curr.index] +
                                    (curr.maha2.2 - maha2.2[curr.index])*
                                    (weights[curr.index+1] - weights[curr.index])/
                                    (maha2.2[curr.index+1]-maha2.2[curr.index])
                              } else {
                                 ## If not, will use 'get.weights' for this maha.
                                 notInterpol[notInterpolcounter] <- ind
                                 notInterpolcounter <- notInterpolcounter + 1
                              }
                           } else {
                              curr.index <- curr.index + 1
                           }
                        }
                     }
                  }
                  ## Now need to approximate weights for those maha in 'notInterpol'
                  if(notInterpolcounter > 1) {
                     notInterpol <- notInterpol[1:(notInterpolcounter-1)]
                     weights.obj <-
                        weights_(maha2.2.new[notInterpol], qW = qW, nu = nu.est,
                                 lrdet = lrdet, d = d, special.mix = special.mix,
                                 control = control, verbose = FALSE)
                     weights.new[notInterpol] <- weights.obj$weights
                     weights.warn.count <-
                        weights.warn.count + weights.obj$weights.warn.count
                     ## Add estimated weights to 'weights' and corresponding
                     ## 'maha2.2.new' to 'maha2.2' so that they can be reused
                     maha2.2 <- c(maha2.2, maha2.2.new[notInterpol])
                     temp.ordering <- order(maha2.2) # only needed here
                     maha2.2 <- maha2.2[temp.ordering] # sort 'maha2.2' again
                     weights <- c(weights, weights.new[notInterpol])[temp.ordering] # and weights accordingly
                     length.maha <- length.maha + notInterpolcounter - 1
                  }
                  ## Recover original ordering and set negative weights to zero.
                  weights.new <- weights.new[order(order.maha2.2.new)] ## TODO: Omit?
                  if(verbose >= 3) cat(".") # print dot after estimation of weights.
               } # done estimating 'weights.new'

               ## Get new 'scale.est': 1/n * sum_{i=1}^n weights_i (x_i-mu)(x_i-mu)^T
               ## where 'mu' corresponds to current 'loc.est'
               scale.est.new <- if(do.scale){
                  crossprod(sqrt(weights.new)*sweep(x, 2, loc.est,
                                                    check.margin = FALSE))/n
               } else scale.est
               ## Get new 'loc.est': sum_{i=1}^n weights_i x_i / (sum weights)
               ## as.vector because we need 'loc.est' as a vector, not (d, 1) matrix
               loc.est.new <- if(do.loc){
                  as.vector(crossprod(x, weights.new)/sum(weights.new))
               } else loc.est
               ## Check convergence
               scale.est.rel.diff <- abs((scale.est - scale.est.new)/scale.est)
               loc.est.rel.diff   <- abs((loc.est - loc.est.new)/loc.est)
               converged.locscale <- if(do.scale & do.loc){
                  (max(loc.est.rel.diff) < control$ECME.rel.conv.tol[2]) &
                     (max(scale.est.rel.diff) < control$ECME.rel.conv.tol[2])
               } else if (do.loc) {
                  (max(loc.est.rel.diff) < control$ECME.rel.conv.tol[2])
               } else if(do.scale){
                  (max(scale.est.rel.diff) < control$ECME.rel.conv.tol[2])
               }
               ## Update counter
               iter.locscaleupdate <- iter.locscaleupdate + 1
               ## Update 'loc.est' and 'scale.est'
               loc.est     <- loc.est.new
               scale.est   <- scale.est.new
            } # done updating 'loc.est' and 'scale.est'
            if(verbose >= 3) cat(paste0(".DONE (", iter.locscaleupdate, " iterations needed)", '\n'))
         }

         ## 2.2  Update 'nu.est' with 'loc.est' and 'scale.est' held fixed #####

         if(control$ECMEstep.do.nu) {
            ## New subsample used for this 'nu' update?
            if(control$resample && size.subsample < n) {
               set.seed(NULL) # destroy potentially resetted seed
               sampled.ind <- sample(n, size.subsample)
               tx.sub      <- tx[,sampled.ind, drop = FALSE]
            }
            ## Optimize neg.log.likelihood over 'nu'
            if(verbose >= 3)
               cat(paste0("    Optimizing likelihood over 'nu' with new 'loc' and 'scale'"))
            est.obj <-
               estim_nu(tx, qW = qW, init.nu = nu.est, loc = loc.est,
                        scale = scale.est, mix.param.bounds = mix.param.bounds,
                        special.mix = special.mix, control = control,
                        control.optim = control$control.optim, verbose = verbose)
            ## Extract results and check convergence
            nu.est.rel.diff <- abs(((nu.est.new <- est.obj$nu.est) - nu.est)/nu.est)
            nu.est <- nu.est.new
            max.ll <- est.obj$max.ll
            dnvmix.warn.count <- dnvmix.warn.count + est.obj$dnvmix.warn.count
            if(verbose >= 3)
               cat(paste0("DONE (", est.obj$ll.counts, " calls to likelihood needed)", '\n'))
         } else {
            nu.est.rel.diff <- 0
         }

         ## 2.3  Check convergence and update various quantities  ##############

         converged <- if(iter.ECME >= control$ECME.miniter) {
            prod(abs(nu.est.rel.diff) < control$ECME.rel.conv.tol[3])
         } else FALSE
         ## Update counter and 'nu.ests.ll'
         iter.ECME <- iter.ECME + 1
         ## Store new 'nu.est' along with log-likelihood
         nu.ests.ll[current.iter.total, ] <- c(nu.est, max.ll)
         ## If 'converged', set all future iteration values to current one
         if(converged) iter.converged <- iter.ECME
         current.iter.total <- current.iter.total + 1

      } # end while()
      ECME.convd <- if(iter.ECME == control$ECME.maxiter & !converged){
         if(verbose >= 1) # print warning
            warning("Maximum number of ECME iterations exhausted, consider increasing 'ECME.maxiter' in the 'control' argument.")
         FALSE
      } else TRUE
   } #end if(control$ECMEstep)

   ## 3 Another last 'nu.est' update with full sample? #########################

   ## This is only relevant/useful when size.subsample < n
   if(control$laststep.do.nu) {
      if(verbose >= 2)
         cat(paste0("Step 3: One last 'nu' update", '\n'))
      if(verbose >= 3)
         cat(paste0("  Optimizing likelihood over 'nu' with new 'loc' and 'scale'"))
      ## One last nu update with the *full* sample and 'control.optim.laststep' as control
      ## (as oppsed to 'control.optim')
      est.obj <- estim_nu(tx, qW = qW, init.nu = nu.est, loc = loc.est,
                          scale = scale.est, mix.param.bounds = mix.param.bounds,
                          special.mix = special.mix, control = control,
                          control.optim = control$control.optim.laststep,
                          verbose = verbose)
      dnvmix.warn.count <- dnvmix.warn.count + est.obj$dnvmix.warn.count
      nu.est <- est.obj$nu.est
      max.ll <- est.obj$max.ll

      ## Store new 'nu.est' along with log-likelihood
      nu.ests.ll[current.iter.total, ] <- c(nu.est, max.ll)
      ## Grab relevant rows of 'nu.ests.ll'
      nu.ests.ll <- nu.ests.ll[1:current.iter.total, , drop = FALSE]
      if(verbose >= 3)
         cat(paste0("DONE (", est.obj$ll.counts, " calls to likelihood needed)", '\n'))
   } else {
      ## Grab relevant rows of 'nu.ests.ll'
      nu.ests.ll <- nu.ests.ll[1:(current.iter.total-1), , drop = FALSE]
   }

   ## 4 Return  ################################################################
   if(verbose >= 2) cat(paste0("RETURN.", '\n'))

   ## Handle warnings
   if(verbose){
      if(dnvmix.warn.count > 0)
         warning(paste0("Error estimation in ", dnvmix.warn.count, " call(s) to the likelihood function was unreliable."))
      if(weights.warn.count > 0)
         warning(paste0("Error estimation in ", weights.warn.count, " update(s) of the weights was unreliable."))
   }

   ## Return S3-class object of type 'fitnvmix'
   return(class_fitnvmix(nu = nu.est, loc = loc.est, scale = scale.est,
                         max.ll = max.ll, data = x, 
                         init.size.subsample = init.size.subsample,
                         size.subsample = size.subsample,
                         dnvmix.warn.count = dnvmix.warn.count,
                         iter.converged = iter.converged, nu.ests.ll = nu.ests.ll,
                         qmix = qmix, is.mvn = is.mvn, is.mvt = is.mvt,
                         do.loc = do.loc, do.scale = do.scale, call = call,
                         ECME.convd = ECME.convd))
}


### S3 class functions and methods #############################################

#' Function to define S3 class 'fitnvmix'
#'
#' @param nu MLE for 'nu'
#' @param scale MLE for 'scale'
#' @param loc MLE for 'loc'
#' @param max.ll maximum log-likelihood at MLEs
#' @param data input data matrix
#' @param init.size.subsample subsample size for initial parameter
#' @param size.subsample subsample size for ECME iterations
#' @param dnvmix.warn.count number of warnings caused by 'dnvmix()'
#' @param weights.warn.count number of warings caused by 'get_weights()'
#' @param iter.converged number of iterations needed until convergence
#' @param nu.ests.ll matrix of estimates of 'nu' with corresponding loglikelihood
#                 in each iteration
#' @param qmix 'qmix' which was passed to 'fitnvmix()'
#' @param is.mvn logical if distribution is multivariate normal
#' @param is.mvt logical if distribution is multivariate t
#' @param do.loc logical if 'loc' is being estimated
#' @param do.scale logical if 'scale' is being estimated
#' @param call language object; function call to 'fitnvmix()'
#' @param ECME.convd logical if ECME converged
#' @return S3 object of class 'fitnvmix'
#' @author Erik Hintz
class_fitnvmix <- function(nu, scale, loc, max.ll, data, init.size.subsample,
                           size.subsample, dnvmix.warn.count = 0, weights.warn.count = 0,
                           iter.converged, nu.ests.ll, qmix, is.mvn = FALSE,
                           is.mvt = FALSE, do.loc = TRUE, do.scale = TRUE,
                           call, ECME.convd){
   res <- if(is.mvn){
      list(nu = NULL, loc = loc, scale = scale, max.ll = max.ll,
           data = data,
           warn.count = list(dnvmix = 0, weights = 0), # no warnings
           iter.converged = 0, nu.ests.ll = NULL, is.mvn = TRUE, is.mvt = FALSE,
           qmix = "constant",
           do.loc = do.loc, do.scale = do.scale,
           init.size.subsample = NULL, size.subsample = NULL,
           call = call, ECME.convd = TRUE)
   } else {
      list(nu = nu, loc = loc, scale = scale, max.ll = max.ll, data = data,
           warn.count = list(dnvmix = dnvmix.warn.count, weights = weights.warn.count),
           iter.converged = iter.converged, nu.ests.ll = nu.ests.ll, is.mvn = FALSE,
           is.mvt = is.mvt, qmix = qmix, do.loc = do.loc, do.scale = do.scale,
           init.size.subsample = init.size.subsample,
           size.subsample = size.subsample, call = call, ECME.convd = ECME.convd)
   }
   ## Return object of class 'fitnvmix'
   structure(res, class = "fitnvmix")
}

## Method 'print' for S3 class 'fitnvmix'
print.fitnvmix <- function(x, ..., digits = max(3, getOption("digits") - 3)){
   ## Grab dimension 'd' and sample size 'n'
   n <- nrow(x$data)
   d <- ncol(x$data)
   ## Print function call to fitnvmix()
   cat("Call: ", deparse(x$call), "\n", sep = "")
   ## Print information about input data
   cat(sprintf(
      "Input data: %d %d-dimensional observations.\n", n, d))
   ## Print information about the distribution (and wether 'loc'/'scale' provided)
   string.provided <- if(!x$do.loc & !x$do.scale){
      "with known 'loc' vector and known 'scale' matrix."
   } else if (!x$do.loc){
      "with known 'loc' vector."
   } else if(!x$do.scale){
      "with known 'scale' matrix."
   } else "with unknown 'loc' vector and unknown 'scale' matrix."
   string.provided.mix <- if(x$is.mvn) "as multivariate normal" else if(x$is.mvt)
      "as multivariate t" else "through quantile function of the mixing variable "
   cat("Normal variance mixture specified", string.provided.mix, "")
   if(!x$is.mvt & !x$is.mvn){
      cat("\n", "   ", deparse(x$qmix), "\n")
   }
   cat(string.provided, "\n")
   ## Print log-likelihood at estimated parameters
   estimated.string <- if(is.character(x$qmix)) "" else "Approximated "
   cat(sprintf("%slog-likelihood at reported parameter estimates: %f \n",
               estimated.string, round(x$max.ll, digits)), sep = "")
   if(!x$is.mvn){
      ## Print subsample and convergence detection
      if(x$size.subsample < n)
         cat(sprintf("Estimation carried out on subsample of size %d",
                     x$size.subsample, ".", "\n"))
      convd.string <- if(x$ECME.convd) "convergence detected." else
         "convergence not detected."
      cat(sprintf("Termination after %d iterations, %s \n", x$iter.converged,
                  convd.string))
      ## Print mixing parameters
      if(x$is.mvt){
         cat("Estimated degrees-of-freedom:", '\n')
         if(any(names(x) == "df")) print(x$df, digits = digits) else
            print(x$nu, digits = digits)
      } else {
         cat("Estimated mixing parameter(s) 'nu':", '\n')
         print(x$nu, digits = digits)
      }
   }
   ## Print estimated 'loc' and 'scale'
   estim.prov.loc <- if(x$do.loc) "Estimated" else "Provided"
   estim.prov.scale <- if(x$do.scale) "Estimated" else "Provided"
   cat(estim.prov.loc, "'loc' vector: ", '\n')
   print(x$loc, digits = digits)
   cat(estim.prov.scale, "'scale' matrix: ", '\n')
   print(x$scale, digits = digits)
   invisible(x) # return
}


## Method 'summary' for S3 class 'fitnvmix'
summary.fitnvmix <- function(object, ..., digits = max(3, getOption("digits") - 3)){
   ## Grab dimension 'd' and sample size 'n'
   n <- nrow(object$data)
   d <- ncol(object$data)
   ## Print function call to fitnvmix()
   cat("Call: ", deparse(object$call), "\n", sep = "")
   ## Print information about input data
   cat(sprintf(
      "Input data: %d %d-dimensional observations.\n", n, d))
   ## Print information about the distribution (and wether 'loc'/'scale' provided)
   string.provided <- if(!object$do.loc & !object$do.scale){
      "with known 'loc' vector and known 'scale' matrix."
   } else if (!object$do.loc){
      "with known 'loc' vector."
   } else if(!object$do.scale){
      "with known 'scale' matrix."
   } else "with unknown 'loc' vector and unknown 'scale' matrix."
   string.provided.mix <- if(object$is.mvn) "as multivariate normal" else if(object$is.mvt)
      "as multivariate t" else "through quantile function of the mixing variable "
   cat("Normal variance mixture specified", string.provided.mix, "")
   if(!object$is.mvt & !object$is.mvn){
      cat("\n", "   ", deparse(object$qmix), "\n")
   }
   cat(string.provided, "\n")
   ## Print log-likelihood at estimated parameters
   estimated.string <- if(is.character(object$qmix)) "" else "Approximated "
   cat(sprintf("%slog-likelihood at reported parameter estimates: %f \n",
               estimated.string, round(object$max.ll, digits)), sep = "")
   if(!object$is.mvn){
      ## Print subsample and convergence detection
      if(object$size.subsample < n)
         cat(sprintf("Estimation carried out on subsample of size %d",
                     object$size.subsample, ".", "\n"))
      convd.string <- if(object$ECME.convd) "convergence detected." else
         "convergence not detected."
      cat(sprintf("Termination after %d iterations, %s \n", object$iter.converged,
                  convd.string))
      ## Print mixing parameters
      if(object$is.mvt){
         cat("Estimated degrees-of-freedom:", '\n')
         if(any(names(object) == "df")) print(object$df, digits = digits) else
            print(object$nu, digits = digits)
      } else {
         cat("Estimated mixing parameter(s) 'nu':", '\n')
         print(object$nu, digits = digits)
      }
   }
   ## Print estimated 'loc' and 'scale'
   estim.prov.loc <- if(object$do.loc) "Estimated" else "Provided"
   estim.prov.scale <- if(object$do.scale) "Estimated" else "Provided"
   cat(estim.prov.loc, "'loc' vector: ", '\n')
   print(object$loc, digits = digits)
   cat(estim.prov.scale, "'scale' matrix: ", '\n')
   print(object$scale, digits = digits)
   
   ## -- up to here same as print.fitnvmix() --
   
   cat("\n")
   if(object$init.size.subsample < n)
      cat("Initial estimate obtained from subsample of size ",
          object$init.size.subsample, ". \n", sep = "")
   print(object$nu.ests.ll)
   ## Print information about warnings
   if(any(object$warn.count > 0)){
      cat("Error estimation in ")
      if(object$warn.count$dnvmix > 0){
         cat(paste0(object$warn.count$dnvmix, " call(s) to the likelihood function "))
         if(object$warn.count$weights > 0)
            cat("and in ")
      }
      if(object$warn.count$weights > 0){
         cat(paste0(object$warn.count$weights, " update(s) of the weights "))
      }
      cat("may have been unreliable. \n")
   }
   invisible(object) # return
}

## Method 'plot' for S3 class 'fitnvmix'
plot.fitnvmix <- function(x, ...){
   if(x$is.mvn){
      cat("Nothing to plot in the case of a multivariate normal distribution.")
   } else {
      ## Grab length of 'nu' and ranges for the two y-axes
      length.par <- if(x$is.mvt){
         ylab = "Estimated df"
         1
      } else {
         ylab <- expression(hat(nu))
         length(x$nu)
      }
      y_rg_nu <- range(x$nu.ests.ll[, 1:length.par]) # range of estimates
      y_rg_ll <- range(x$nu.ests.ll[, length.par + 1]) # range of log-likelihood
      iters <- 1:nrow(x$nu.ests.ll) # iterations
      ## Prepare legend
      lgnd <- if(x$is.mvt) c("log-likelihood", "df") else {
         tmp <- vector("expression", length.par + 1)
         tmp[[1]] <- "log-likelihood"
         for(i in 1:length.par)
            tmp[[i+1]] <- bquote(hat(nu)[.(i)])
         tmp
      }
      ## Plot
      def.par <- par(no.readonly = TRUE) # save default, for resetting...
      par(mar = c(4, 3, 3, 3) + 0.15)
      ## Plot estimates as a function of iterations
      plot(NA, xlab = "Iteration", ylab = "", ylim = y_rg_nu,
           xlim = range(iters-1), xaxt = "n", yaxt = "n")
      for(i in 1:length.par)
         lines(iters-1, x$nu.ests.ll[, i], col = i+1, lty = i+1)
      axis(2, ylim = y_rg_nu, lwd = 1)
      mtext(2, text = ylab, line = 1.9)
      ## Plot log-likelihood, axes and legend
      par(new = T)
      plot(iters-1, c(NA, x$nu.ests.ll[-1, length.par + 1]), type = 'l',
           axes = F, xlab = "", ylab = "")
      axis(4, ylim = y_rg_ll, lwd = 1, line = 0)
      mtext(4, text = "log-likelihood", line = 2)
      axis(1, pretty(range(iters-1), length(iters)))
      legend("topright", legend = lgnd, lty = 1:(length.par + 1),
             col = 1:(length.par+1), bty = 'n')
      par(def.par) # reset to default
   }
   invisible(x)
}
