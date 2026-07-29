### pgammamix() ################################################################

##' @title Distribution Function of the Mahalanobis Distance of a Normal
##'        Variance Mixture
##' @param x n-vector of evaluation points
##' @param qmix see ?pnvmix()
##' @param lower.tail logical if P(X <= m) or P(X>m) shd be returned
##' @param d dimension of the underlying normal variance mixture
##' @param control see ?pnvmix()
##' @param verbose see ?pnvmix()
##' @param ... asee ?pnvmix()
##' @return n-vector with computed cdf values and attributes 'error'
##'         (error estimate) and 'numiter' (number of while-loop iterations)
##' @author Erik Hintz and Marius Hofert
pgammamix <- function(x, qmix, d, lower.tail = TRUE,
                      control = list(), verbose = TRUE, ...)
{
    ## Checks
    if(!is.vector(x)) x <- as.vector(x)
    n <- length(x) # length of input
    ## Deal with algorithm parameters, see also get_set_param():
    ## get_set_param() also does argument checking, so not needed here.
    control <- get_set_param(control)

    ## 1 Define the quantile function of the mixing variable ###################
    mix_list      <- get_mix_(qmix = qmix, callingfun = "pgammamix", ... ) 
    qW            <- mix_list[[1]] # function(u)
    special.mix   <- mix_list[[2]]
    
    ## Build result object
    pres <- rep(0, n) # n-vector of results
    abserror <- rep(0, n)
    relerror <- rep(0, n) 
    notNA <- which(!is.na(x))
    pres[!notNA] <- NA
    x <- x[notNA, drop = FALSE] # non-missing data (rows)
    ## Counter
    numiter <- 0 # initialize counter (0 for 'inv.gam' and 'is.const.mix')
    ## Deal with special distributions
    if(!is.na(special.mix)) {
        if(!(special.mix == "pareto")) {
            ## Only for "inverse.gamma" and "constant" do we have analytical forms:
            pres[notNA] <- switch(special.mix,
                                  "inverse.gamma" = {
                                      ## D^2 ~ d* F(d, nu)
                                      pf(x/d, df1 = d, df2 = mix_list$param,
                                         lower.tail = lower.tail)
                                  },
                                  "constant" = {
                                      ## D^2 ~ chi^2_d = Gamma(shape = d/2, scale = 2)
                                      pgamma(x, shape = d/2, scale = 2,
                                             lower.tail = lower.tail)
                                  })
            ## Return in those cases
            attr(pres, "abs. error") <- rep(0, n)
            attr(pres, "rel. error") <- rep(0, n)
            attr(pres, "numiter") <- numiter
            return(pres)
        }
    }
    ## Define various quantites for the RQMC procedure
    dblng           <- (control$increment == "doubling")
    B               <- control$B # number of randomizations
    current.n       <- control$fun.eval[1] #initial sample size
    total.fun.evals <- 0
    ZERO            <- .Machine$double.neg.eps
    ## Absolte/relative precision?
    if(is.na(control$pgammamix.reltol)) {
        ## Use absolute error
        tol <- control$pgammamix.abstol
        do.reltol <- FALSE
    } else {
        ## Use relative error
        tol <- control$pgammamix.reltol
        do.reltol <- TRUE
    }
    ## Store seed if 'sobol' is used to get the same shifts later
    if(control$method == "sobol") {
       seeds_ <- sample(1:(1e3*B), B) # B seeds for 'sobol()'
    }
    ## Additional variables needed if the increment chosen is "dblng"
    if(dblng) {
        if(control$method == "sobol") useskip <- 0
        denom <- 1
    }
    ## Matrix to store RQMC estimates
    rqmc.estimates <- matrix(0, ncol = n, nrow = B)
    CI.factor.sqrt.B <- control$CI.factor / sqrt(B) # needed again and again
    ## Initialize 'max.error' to > tol so that we can enter the while loop
    max.error <- tol + 42

    ## 2 Main loop #############################################################

    ## while() runs until precision abstol is reached or the number of function
    ## evaluations exceed fun.eval[2]. In each iteration, B RQMC estimates of
    ## the desired probability are calculated.
    while(max.error > tol && numiter < control$max.iter.rqmc &&
          total.fun.evals < control$fun.eval[2])
    {

        for(b in 1:B) {

            ## 2.1 Get the point set
            U <- switch(control$method,
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
                        }) # sorted for later!

            ## 2.2 Evaluate the integrand at the (next) point set
            W <- pmax(qW(U), ZERO) # realizations of the mixing variable
            next.estimate <- .colMeans(pgamma(
                sapply(x, function(i) i/W), shape = d/2, scale = 2,
                lower.tail = lower.tail), current.n, n, 0)

            ## 2.3 Update RQMC estimates
            rqmc.estimates[b,] <-
                if(dblng) {
                    ## In this case both, rqmc.estimates[b,] and
                    ## next.estimate depend on n.current points
                    (rqmc.estimates[b,] + next.estimate)/denom
                } else {
                    ## In this case, rqmc.estimates[b,] depends on
                    ## numiter * n.current points whereas next.estimate
                    ## depends on n.current points
                    (numiter * rqmc.estimates[b,] + next.estimate)/(numiter + 1)
                }
        } # end for(b in 1:B)

        ## Update of various variables
        ## Double sample size and adjust denominator in averaging as well as useskip
        if(dblng) {
            ## Change denom and useksip (exactly once, in the first iteration)
            if(numiter == 0) {
                denom <- 2
                useskip <- 1
            } else {
                ## Increase sample size n. This is done in all iterations
                ## except for the first two
                current.n <- 2 * current.n
            }
        }
       ## Total number of function evaluations
       total.fun.evals <- total.fun.evals + B * current.n
       numiter <- numiter + 1
       ## Update error. The following is slightly faster than 'apply(..., 2, var)'
       pres[notNA] <- .colMeans(rqmc.estimates, B, n , 0)
       vars <- .colMeans((rqmc.estimates - rep(pres, each = B))^2, B, n, 0)
       error <- if(!do.reltol) {
          sqrt(vars)*CI.factor.sqrt.B
       } else {
          sqrt(vars)/pres*CI.factor.sqrt.B
       }
       max.error <- max(error)
    } # while()
    if(verbose & max.error > tol)
       warning("Tolerance not reached for all inputs; consider increasing 'max.iter.rqmc' in the 'control' argument.")
    
    ## Compute absolute and relative errors (once here at the end)
    abserror[notNA] <- if(do.reltol){
       relerror[notNA] <- error
       relerror[notNA] * pres[notNA] 
    } else { # error is absolute error
       relerror[notNA] <- error / pres[notNA] 
       error 
    }
    ## Return
    attr(pres, "abs. error") <- abserror
    attr(pres, "rel. error") <- relerror
    attr(pres, "numiter") <- numiter
    pres
}


### qgammamix() ################################################################

##' @title Quantile function of the mahalanobis distance of a normal variance mixture
##' @param u see ?qnvmix()
##' @param qmix see ?qnvmix()
##' @param d dimension of the underlying normal variance mixture
##' @param control see ?qnvmix()
##' @param verbose see ?qnvmix()
##' @param q.only see ?qnvmix()
##' @param stored.values see ?qnvmix()
##' @param ... see ?qnvmix()
##' @return see ?qnvmix()
qgammamix <- function(u, qmix, d, control = list(), verbose = TRUE, q.only = TRUE,
                      stored.values = NULL, ...)
    quantile_(u, qmix = qmix, which = "maha2", d = d, control = control,
                      verbose = verbose, q.only = q.only,
                      stored.values = stored.values, ...)


### rgammamix() ################################################################

##' @title Random Number Generator for Gamma mixtures (= Mahalanobis distances for
##'        normal mixture models)
##' @param n see ?rnvmix()
##' @param rmix ?rnvmix()
##' @param qmix ?rnvmix()
##' @param d dimension of the underlying normal variance mixture, see below under 'return'
##' @param method ?rnvmix()
##' @param skip ?rnvmix()
##' @param ... n-vector or realizations
##' @return n-vector of samples of W*X where W~qmix/rmix and X~chi^2_d = Gamma(shape = d/2, scale = 2)
##' @author Marius Hofert and Erik Hintz
rgammamix <- function(n, rmix, qmix, d, method = c("PRNG", "sobol", "ghalton"), 
                      skip = 0, ...)
    rnvmix_(n, rmix = rmix, qmix = qmix, loc = rep(0, d), scale = diag(d),
            factor = diag(d), method = method, skip = skip, which = "maha2", ...)


### dgammamix() ################################################################

##' @title Density function of the mahalanobis distance of a normal variance mixture
##' @param x n-vector of evaluation points
##' @param qmix see ?pnvmix()
##' @param d dimension of the underlying normal variance mixture
##' @param control see ?pnvmix()
##' @param verbose see ?pnvmix()
##' @param log logical if log-density shall be returned
##' @param ... see ?pnvmix()
##' @return n-vector with computed density values and attributes 'error'
##'         (error estimate) and 'numiter' (number of while-loop iterations)
##' @author Erik Hintz and Marius Hofert
dgammamix <- function(x, qmix, d, control = list(), verbose = TRUE, log = FALSE, ...)
{
    ## Checks
    if(!is.vector(x)) x <- as.vector(x)
    n <- length(x) # dimension
    stopifnot(d >= 1)
    verbose <- as.logical(verbose)
    ## Deal with algorithm parameters, see also get_set_param()
    control <- get_set_param(control)
    ## Build result object (log-density)
    lres  <- rep(-Inf, n) # n-vector of results
    abserror <- rep(NA, n)
    relerror <- rep(NA, n) 
    notNA <- which(!is.na(x))
    lres[!notNA] <- NA
    x <- x[notNA, drop = FALSE] # non-missing data (rows)
    ## Define the quantile function of the mixing variable 
    if(!hasArg(qmix)) qmix <- NULL # needed for 'get_mix_()'
    mix_list      <- get_mix_(qmix = qmix, callingfun = "dgammamix", ... ) 
    qW            <- mix_list[[1]] # function(u)
    special.mix   <- mix_list[[2]] # NA or string 
    ## Counter
    numiter <- 0 # initialize counter (0 for 'inv.gam' and 'is.const.mix')
    ## Deal with special distributions
    if(!is.na(special.mix)) {
        if(!(special.mix == "pareto")) {
            ## Only for "inverse.gamma" and "constant" do we have analytical forms
            lres[notNA] <- switch(special.mix,
                                  "inverse.gamma" = {
                                      ## D^2 ~ d* F(d, nu)
                                      df(x/d, df1 = d, df2 = mix_list$param, 
                                         log = TRUE) - log(d)
                                  },
                                  "constant" = {
                                      ## D^2 ~ chi^2_d = Gamma(shape = d/2, scale = 2)
                                      dgamma(x, shape = d/2, scale = 2, log = TRUE)
                                  })
            ## Return in those cases
            attr(lres, "abs. error") <- rep(0, n)
            attr(lres, "rel. error") <- rep(0, n)
            attr(lres, "numiter") <- numiter
            if(log) return(lres) else return(exp(lres))
        }
    }
    ## General case of a multivariate normal variance mixture (=> RQMC procedure)
    ## Prepare inputs for 'densmix_()'
    ## Sort input 'x' increasingly and store ordering for later
    ordering.x <- order(x)
    x.ordered  <- x[ordering.x]
    ## Define log-constant for the integration
    lconst <- -lgamma(d/2) - d/2 * log(2) + (d/2-1)*log(x.ordered)
    ## Call 'densmix_()' (which itself calls C-Code and handles warnings):
    est.list <- densmix_(qW, maha2.2 = x.ordered/2, lconst = lconst,
                             d = d, control = control, verbose = verbose)
    ## Grab results, correct 'error' and 'lres' if 'log = FALSE'
    ldens <- est.list$ldensities[order(ordering.x)]
    abserror[notNA] <- est.list$abserror[order(ordering.x)]
    relerror[notNA] <- est.list$relerror[order(ordering.x)]
    ## Correct results and error if 'log = FALSE'
    if(!log){
       ldens <- exp(ldens)
       ## CI for mu: exp(logmu_hat +/- abserr(logmu_hat))) = (lower, upper)
       ## => compute max. error on mu_hat as max( (upper - mu), (mu - lower) ) 
       relerror[notNA] <- max( (exp(abserror[notNA]) - 1), (1 - exp(-abserror[notNA])) )
       abserror[notNA] <- ldens * relerror[notNA] # ldens already exponentiated 
    }
    lres[notNA] <- ldens 
    
    ## Return
    ## Note that 'lres' was exponentiated already if necessary.
    attr(lres, "abs. error") <- abserror
    attr(lres, "rel. error") <- relerror
    attr(lres, "numiter") <- est.list$numiter
    lres
}
