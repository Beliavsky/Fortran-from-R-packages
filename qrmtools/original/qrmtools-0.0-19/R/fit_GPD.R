### Fitting GPD(shape, scale) ##################################################

## ## Other packages
## - QRM:
##   + fit.GPD() has various options but, per default, uses probability-weighted
##     moments to get initial values and then MLE based on optim() *with*
##     provided gradient and observed information
## - Renext:
##   + fGPD(); uses a bit of a weird procedure based on the special cases
##     lomax and maxlo
##   + detailed in "Renext Computing Details" (but not available)
## - evd:
##   + fpot() -> fpot.norm()
##   + based on optim() with start values shape = 0, scale = 0
## - evir:
##   + ./R -> pot.R -> gpd(): based on optim() with initial values
##     s2 <- var(excess)
##     shape0 <- -0.5 * (((xbar * xbar)/s2) - 1)
##     scale0 <- 0.5 * xbar * (((xbar * xbar)/s2) + 1)
##     theta <- c(shape0, scale0)
##   + with hessian and observed information
## - fExtremes:
##   + gpdFit() -> .gpdmleFit(): as evir
## - ismev:
##   + gpd.fit(): uses optim(), hardcoded initial values
##     (shape = 0.1, scale = sqrt(6 * var(xdat))/pi)
## - texmex:
##   + evm() -> evm.default() -> evmFit(): based on optim(); uses
##     'family$start(data)' if start not provided (see also .pdf)
##   + 'start()' is a function and for the GPD found in ./R/gpd.R
##   + start() returns a longer vector (unclear why) but seems to use log(mean())
##     as initial value for shape and 1e-05 for scale; quite unclear


### Method-of-moments estimator ################################################

##' @title Method-Of-Moments Estimator
##' @param x numeric vector of data. In the peaks-over-threshold method,
##'        these are the excesses over a sufficiently high threshold.
##' @return numeric(2) with estimates of shape and scale
##' @author Marius Hofert
##' @note See Hosking and Wallis (1987) and Landwehr, Matalas, Wallis (1979)
fit_GPD_MOM <- function(x)
{
    loc.hat <- mean(x)
    sig2.hat <- var(x)
    shape.hat <- (1-loc.hat^2/sig2.hat)/2
    c(shape = shape.hat, scale = max(loc.hat*(1-shape.hat), .Machine$double.eps)) # shape, scale (> 0)
}


### Probability weighted moments estimator #####################################

##' @title Probability Weighted Moments Estimator
##' @param x numeric vector of data. In the peaks-over-threshold method,
##'        these are the excesses over a sufficiently high threshold.
##' @return numeric(2) with estimates of shape and scale
##' @author Marius Hofert
##' @note See Hosking and Wallis (1987) and Landwehr, Matalas, Wallis (1979)
##'       Note that their '-k' and 'alpha' are our 'xi' and 'beta'.
fit_GPD_PWM <- function(x)
{
    ## a_s = estimator (unbiased according to Landwehr, Matalas, Wallis (1979)) / a sample version of M_{1,0,s} = E(X (1-F(X))^s)
    x. <- sort(as.numeric(x)) # CAUTION: If is.xts(x), then sort won't do anything!
    n <- length(x.)
    k <- 1:n
    a0 <- mean(x)
    a1 <- mean(x. * (n-k)/(n-1))
    ## Alternatively (see QRM and Hosking and Wallis (1987)): a1 <- mean(x. * (1-p)) where p <- (k-0.35)/n
    ## Estimators of shape and scale based on a0 and a1
    c(shape = 2 - a0/(a0 - 2 * a1), scale = max((2 * a0 * a1) / (a0 - 2 * a1), .Machine$double.eps)) # shape, scale (> 0)
}


### Maximum likelihood estimator ###############################################

##' @title Log-likelihood of the GPD
##' @param param numeric(2) giving shape and scale (all real here;
##'        if scale <= 0 dGPD(, log = TRUE) returns -Inf)
##' @param x numeric vector of data. In the peaks-over-threshold method,
##'        these are the excesses over a sufficiently high threshold.
##' @return log-likelihood at shape and scale
##' @author Marius Hofert
logLik_GPD <- function(param, x)
    sum(dGPD(x, shape = param[1], scale = param[2], log = TRUE))

##' @title MLE for GPD Parameters
##' @param x numeric vector of data. In the peaks-over-threshold method,
##'        these are the excesses over a sufficiently high threshold.
##' @param init numeric(2) giving the initial values for shape and scale.
##'        Findings of Hosking and Wallis (1987):
##'        - MLE requires n ~>= 500 to be efficient
##'        - MoM reliable for shape > 0.2
##'        - PWM recommended for shape > 0.
##' @param estimate.cov logical indicating whether the asymptotic covariance
##'        matrix of the parameter estimators is to be estimated
##'        (inverse of observed Fisher information (negative Hessian
##'        of log-likelihood evaluated at MLE))
##' @param control see ?optim
##' @param ... additional arguments passed to the underlying optim()
##' @return list with the return value of optim() with the estimated asymptotic
##'         covariance matrix of the parameter estimators
##'         (Cramer--Rao bound = inverse of Fisher information = negative Hessian)
##'         of the parameter estimators appended if estimate.cov.
##' @author Marius Hofert
##' @note - Note that for no initial shape is the support of the density IR and thus
##'         the log-likelihood could be -Inf if (some) x are out of the support.
##'         We avoid this here by a) requiring x >= 0 and b) adjusting scale so
##'         that the log-likelihood is finite at the initial value.
##'       - Most of the results generalize to the case where the GPD has another
##'         parameter for the location (like mu for the GEV distribution). If this
##'         mu would be chosen <= min(x), then one could probably guarantee a
##'         finite log-likelihood.
##'       - Similar to copula:::fitCopula.ml()
##'       - *Expected* information is available, too; see QRM::fit.GPD()
fit_GPD_MLE <- function(x, init = c("PWM", "MOM", "shape0"), estimate.cov = TRUE, control = list(), ...)
{
    ## Checks
    isnum <- is.numeric(init)
    stopifnot(is.numeric(x), is.character(init) || (isnum && length(init) == 2),
              is.logical(estimate.cov), is.list(control))
    if(any(x < 0))
        stop("The support for the two-parameter GPD distribution is >= 0.")

    ## Initial values
    if(!isnum) {
        switch(match.arg(init),
        "PWM" = {
            init <- fit_GPD_PWM(x)
            ## Force initial values to produce finite log-likelihood
            ## (all x need to be in [0, -scale/shape) => choose scale = -shape * max(x) * 1.01)
            if(init[1] < 0) { # shape < 0
                mx <- max(x)
                if(mx >= -init[2]/init[1]) init[2] <- -init[1] * mx * 1.01 # scale = -shape * max(x) * 1.01
            }
        },
        "MOM" = {
            init <- fit_GPD_MOM(x)
            ## Force initial values to produce finite log-likelihood
            ## (all x need to be in [0, -scale/shape) => choose scale = -shape * max(x) * 1.01)
            if(init[1] < 0) { # shape < 0
                mx <- max(x)
                if(mx >= -init[2]/init[1]) init[2] <- -init[1] * mx * 1.01 # scale = -shape * max(x) * 1.01
            }
        },
        "shape0" = {
            ## Idea: Use shape = 0
            init <- c(0, mean(x)) # mean for shape = 0 is scale
        },
        stop("Wrong 'init'"))
    }

    ## Fit
    control <- c(as.list(control), fnscale = -1) # maximization (overwrites possible additionally passed 'fnscale')
    fit <- optim(init, fn = function(param) logLik_GPD(param, x = x),
                 hessian = estimate.cov, control = control, ...)
    names(fit$par) <- c("shape", "scale")
    ## Note: Could incorporate the gradient of the log-likelihood (w.r.t. shape, scale)
    ##       for the methods "BFGS", "CG" and "L-BFGS-B".
    ##       It is: (0, (-n+sum(x)/scale)/scale) for shape != 0 and
    ##       ( (1/shape^2)*sum(log1p(shape*x/scale))-(shape+1)*sum(shape*x/(scale+shape*x)),
    ##         -n/scale + (1+1/shape) * sum((shape*x/scale)/(scale+shape*x)) )

    ## Estimate of the asymptotic covariance matrix and standard errors
    ## of the parameter estimators
    ## Note: manually:
    ##       fisher <- -hessian(logLik_GPD(fit$par, x = x)) # observed Fisher information (estimate of Fisher information)
    ##       Cov <- solve(fisher)
    ##       std.err <- sqrt(diag(Cov)) # see also ?fit_GPD_MLE
    if(estimate.cov) {
        negHessianInv <- catch(solve(-fit$hessian))
        if(is(negHessianInv$error, "error")) {
            warning("Hessian matrix not invertible: ", negHessianInv$error)
            Cov <- matrix(NA_real_, 0, 0)
            SE <- numeric(0)
        } else {
            Cov <- negHessianInv$value # result on warning or on 'worked'
            rownames(Cov) <- c("shape", "scale")
            colnames(Cov) <- c("shape", "scale")
            SE <- sqrt(diag(Cov))
            names(SE) <- c("shape", "scale")
        }
    } else {
        Cov <- matrix(NA_real_, 0, 0)
        SE <- numeric(0)
    }

    ## Return (could create an object here)
    c(fit, list(Cov = Cov, SE = SE))
}
