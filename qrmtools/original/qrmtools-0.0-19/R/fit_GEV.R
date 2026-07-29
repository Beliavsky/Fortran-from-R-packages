### Fitting GEV(shape, loc, scale) #############################################

## ## Other packages
## - QRM:
##   + fit.GEV() uses hard-coded value for shape (= 0.1) and loc and scale
##     from case shape = 0 (see below)
## - Renext:
##   + also based on optim()
##   + fGEV.MAX() -> calls parIni.MAX(); detailed in "Renext Computing Details"
##     (but not available)
##   + some sort of weird procedure based on a regression idea of some sort
##   + good 'caution': log-lik = Inf for shape <= -1 but optimization
##     done in unconstrained way. Typically shape > -1, but if shape <= -1, then
##     meaningless.
##     For shape < -0.5, treat with care.
## - evd:
##   + also based on optim()
##   + fextreme(): requires 'start' to be provided
##   + fgev(): requires 'start' to be provided
## - evir:
##   + ./R -> bmax.R -> gev(): hardcoded shape = 0.1 as 'QRM'
## - extRemes:
##   + fevd(): huge, uses optim() with 'init.pars', can provide 'initial' or
##     determine them
##   + 'initial' is determine (on 'find.init') via L moments or moments and
##     hard-coded (shape = 0.01, loc = 0, sig = 1) on failure of both methods
## - fExtremes:
##   + gevFit() -> .gevFit() -> .gevmleFit(): uses optim(), hardcoded initial
##     values as 'QRM' (refers to evir for that)
## - ismev:
##   + gev.fit(): uses optim(), hardcoded initial values as in 'QRM'
## - lmom:
##   + pelp(): based on optim(), 'start' needs to be provided
## - texmex:
##   + evm() -> evm.default() -> evmFit(): based on optim(); uses
##     'family$start(data)' if start not provided (see also .pdf)
##   + 'start()' is a function and for the GEV found in ./R/gev.R
##   + start() returns a longer vector (unclear why) but seems to use mean()
##     as initial value for loc and log(IQR(data)/2) as initial value for shape;
##     quite unclear
##   + points out that MLE will often fail with small sample size


### Quantile based estimator ###################################################

##' @title Quantile Based Estimator
##' @param x numeric vector of data. In the block maxima method, these are the
##'        block maxima (based on block size n).
##' @param p numeric(3) specifying the probabilities whose quantiles are
##'        matched
##' @param cutoff value after which exp(-x) is truncated to 0
##' @return numeric(3) with estimates of shape, loc, scale
##' @author Marius Hofert
##' @note Determines shape, loc, scale such that the empirical p-quantiles are matched
fit_GEV_quantile <- function(x, p = c(0.25, 0.5, 0.75), cutoff = 3)
{
    stopifnot(length(p) == 3, 0 < p, p < 1, cutoff > 0)
    q <- quantile(x, probs = p, names = FALSE) # empirical p-quantiles
    if(length(unique(q)) < 3)
        stop("quantile(x, probs = p) does not return unique quantiles.")
    q.diff <- diff(q)
    y <- q.diff[1]/q.diff[2] # (q[2] - q[1]) / (q[3] - q[2]) = (H^{-1}(p2) - H^{-1}(p1)) / (H^{-1}(p3) - H^{-1}(p2))
    l <- log(-log(p)) # decreasing in p
    a <- rev(diff(rev(l))) # l[1] - l[2]; l[2] - l[3]
    a. <- a[1]/a[2]
    ## Initial value for shape
    shape.hat <- if(y < 1/expm1(cutoff/a.)) {
                  log1p(1/y)/a[2]
              } else if(y <= a.) {
                  m1 <- a[2]/cutoff * log(a./(expm1(a.*cutoff)))
                  log(y/a.) / m1
              } else if(y <= expm1(a.*cutoff)) {
                  m2 <- -a[1]/cutoff * log(a.*expm1(cutoff/a.))
                  log(y/a.) / m2
              } else -log1p(y)/a[1]
    ## Initial value for scale
    scale.hat <- if(shape.hat == 0) {
                   q.diff[1] / (-l[2] + l[1])
               } else {
                   shape.hat * q.diff[1] / ((-log(p[2]))^(-shape.hat) - (-log(p[1]))^(-shape.hat))
               }
    ## Initial value for loc
    loc.hat <- if(shape.hat == 0) {
                  q[2] + scale.hat * l[2]
              } else {
                  q[2] - scale.hat/shape.hat * ((-log(p[2]))^(-shape.hat)-1)
              }
    ## Return
    c(shape = shape.hat, loc = loc.hat, scale = scale.hat)
}


### Probability weighted moments estimator #####################################

##' @title Probability Weighted Moments Estimator
##' @param x numeric vector of data. In the block maxima method, these are the
##'        block maxima (based on block size n).
##' @return numeric(3) with estimates of shape, loc, scale
##' @author Marius Hofert
##' @note See Hosking, Wallis, Wood (1985) and Landwehr and Wallis (1979)
##'       Note that their '-k', 'xi' and 'alpha' are our 'xi', 'mu' and 'sigma'.
fit_GEV_PWM <- function(x)
{
    ## b_r = estimator (unbiased according to Landwehr and Wallis (1979)) / a sample version of M_{1,r,0} = E(X F(X)^r)
    x. <- sort(as.numeric(x)) # CAUTION: If is.xts(x), then sort won't do anything!
    n <- length(x.)
    k <- 1:n
    b0 <- mean(x)
    b1 <- mean(x. * (k-1)/(n-1))
    b2 <- mean(x. * (k-1)*(k-2)/((n-1)*(n-2)))
    ## Match probability weighted moments
    ## 1) Theoretically:
    ## shape <- seq(-10, 10, length.out = 65)
    ## y <- (3^shape-1)/(2^shape-1)
    ## plot(shape, y, type = "l") # => starts at 1 (shape = -Inf), increases exponentially)
    y <- (3*b2-b0) / (2*b1-b0) # evaluation point of h^{-1} for h(shape) = (3^shape-1)/(2^shape-1) (>= 1)
    if(y <= 1)
        stop("Cannot invert auxiliary function h involved. Use a different estimator.")
    ## shape.hat <- uniroot(function(x) (3^x-1)/(2^x-1) - y, interval = c(-100, 100))$root
    ## 2) Hosking, Wallis, Wood (1985, (14)) suggest to approximate h^{-1}(y) via g^{-1}(c)
    ##    for c = (2*b1-b0)/(3*b2-b0) - log(2)/log(3) and g^{-1}(c) = -(7.8590 + 2.9554*c)*c.
    ##    This is equivalent to h^{-1}(y) = a/y^2 + b/y + c for a, b, c as follows.
    a <- -2.9554
    b <- -4.1297 # = -(7.8590-2.9554*log(4)/log(3))
    c <- 3.782014 # = -(-7.8590*log(2)/log(3)+2.9554*(log(2)/log(3))^2)
    shape.hat <- a/y^2 + b/y + c
    ## Now compute the estimators for loc and scale
    scale.hat <- if(shape.hat == 0) {
                   (2 * b1 - b0) / log(2)
               } else {
                   (2 * b1 - b0) * shape.hat / (gamma(1 - shape.hat) * (2^shape.hat - 1)) # scale, see Hosking, Wallis, Wood (1985, (15))
               }
    loc.hat <- if(shape.hat == 0) {
                  b0 - scale.hat * 0.5772157 # ... Euler--Mascheroni constant
              } else {
                  b0 - scale.hat * (gamma(1-shape.hat) - 1) / shape.hat # loc, see Hosking, Wallis, Wood (1985, (15))
              }
    ## Return
    c(shape = shape.hat, loc = loc.hat, scale = scale.hat)
}


### Maximum likelihood estimator ###############################################

##' @title Log-likelihood of the GEV
##' @param param numeric(3) giving shape, loc, scale (all real here; if scale <= 0
##'        dGEV(, log = TRUE) returns -Inf)
##' @param x numeric vector of data. In the block maxima method, these are the
##'        block maxima (based on block size n)
##' @return log-likelihood at shape, loc, scale
##' @author Marius Hofert
logLik_GEV <- function(param, x)
    sum(dGEV(x, shape = param[1], loc = param[2], scale = param[3], log = TRUE))

##' @title MLE for GEV Parameters
##' @param x numeric vector of data. In the block maxima method, these are the
##'        block maxima (based on block size n)
##' @param init numeric(3) or string giving the initial values for shape, loc, scale.
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
##' @note 1) similar to copula:::fitCopula.ml()
##'       2) careful for shape <= -0.5 (very short, bounded upper tail):
##'          MLE doesn't have standard asymptotic properties.
##'       3) No method except 'shape0' guarantees that 1 + shape (x - loc) / scale > 0
##'          and thus that the log-likelihood > -Inf for all data points; if
##'          log-likelihood == 0, optim() stops with error:
##'          "Error in optim(init, fn = function(param) logLik_GEV(param, x = x),
##'          hessian = estimate.cov: function cannot be evaluated at initial
##'          parameters"
##'          => below we guarantee to get a finite log-likelihood by doubling scale
##'       4) Other options as default method are 'shape0' with shape = 0 (exact) but
##'          method = "BFGS", or 'PWM"
fit_GEV_MLE <- function(x, init = c("shape0", "PWM", "quantile"),
                        estimate.cov = TRUE, control = list(), ...)
{
    ## Checks
    isnum <- is.numeric(init)
    stopifnot(is.numeric(x), is.character(init) || (isnum && length(init) == 3),
              is.logical(estimate.cov), is.list(control))

    ## Initial values
    if(!isnum) {
        switch(match.arg(init),
        "shape0" = {
            ## Idea: - Use shape = 0 here => loc and scale explicit in terms of mean and
            ##         variance (method-of-moment estimator for loc, scale).
            ##         This guarantees that 1 + shape (x - loc) / scale > 0, so a finite
            ##         log-likelihood.
            ##       - As in 'QRM', 'evir', 'fExtremes', 'ismev'
            ##       - Note that shape = 0 can fail for method = "Nelder-Mead"
            ##         (seen for the Black Monday example)
            ##         => .Machine$double.eps works (but is still unstable, seen for
            ##            GEV_shape_plot() of Black Monday example)
            ##         => use .Machine$double.eps^0.25 as in uniroot()
            scale.hat <- sqrt(6 * var(x)) / pi # var for shape = 0 is (scale \pi)^2 / 6 => scale
            init <- c(.Machine$double.eps^0.25, mean(x) - 0.5772157 * scale.hat, scale.hat) # mean for shape = 0 is loc + sig * gamma => loc; gamma = -digamma(1) = Euler--Mascheroni constant
        },
        "PWM" = {
            init <- fit_GEV_PWM(x)
            ## Force initial values to produce finite log-likelihood
            ## Formerly:
            ## if(method != "shape0") { # only for 'shape0', a non-zero density is guaranteed
            ##     scale.bound <- max(-init[1] * (x - init[2]))
            ##     if(init[3] <= scale.bound)
            ##         init[3] <- 1.05 * scale.bound # blow up by 5%
            ## }
            while(!is.finite(logLik_GEV(init, x = x)) && is.finite(init[3])) init[3] <- init[3] * 2
            ## Note: if !is.finite(init[3]), there's nothing we can do...
        },
        "quantile" = {
            init <- fit_GEV_quantile(x)
            while(!is.finite(logLik_GEV(init, x = x)) && is.finite(init[3])) init[3] <- init[3] * 2
        },
        stop("Wrong 'init'"))
    }

    ## Fit
    control <- c(as.list(control), fnscale = -1) # maximization (overwrites possible additionally passed 'fnscale')
    fit <- optim(init, fn = function(param) logLik_GEV(param, x = x),
                 hessian = estimate.cov, control = control, ...)
    names(fit$par) <- c("shape", "loc", "scale")

    ## Estimate of the asymptotic covariance matrix and standard errors
    ## of the parameter estimators
    ## Note: manually; see QRM::fit.GEV and also copula:::fitCopula.ml():
    ##       fisher <- -hessian(logLik_GEV(fit$par, x = x)) # observed Fisher information (estimate of Fisher information)
    ##       Cov <- solve(fisher)
    ##       std.err <- sqrt(diag(Cov)) # see also ?fit_GEV_MLE
    if(estimate.cov) {
        negHessianInv <- catch(solve(-fit$hessian))
        if(is(negHessianInv, "error")) {
            warning("Hessian matrix not invertible: ", negHessianInv$error)
            Cov <- matrix(NA_real_, 0, 0)
            SE <- numeric(0)
        } else {
            Cov <- negHessianInv$value # result on warning or on 'worked'
            rownames(Cov) <- c("shape", "loc", "scale")
            colnames(Cov) <- c("shape", "loc", "scale")
            SE <- sqrt(diag(Cov))
            names(SE) <- c("shape", "loc", "scale")
        }
    } else {
        Cov <- matrix(NA_real_, 0, 0)
        SE <- numeric(0)
    }

    ## Return (could create an object here)
    c(fit, list(Cov = Cov, SE = SE))
}



