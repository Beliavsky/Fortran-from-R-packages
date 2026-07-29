### d/p/rnvmixcopula() ############################################################

##' @title Density Function of a Multivariate Normal Variance Mixture Copula
##' @param u (n,d) matrix of evaluation points. Have to be in (0,1)
##' @param qmix see ?pnvmix
##' @param scale (d, d)-covariance matrix (scale matrix).
##' @param factor Cholesky factor (lower triangular matrix) of 'scale'
##' @param control see ?get_set_param()
##' @param verbose logical (or integer: 0 = FALSE, 1 = TRUE, 2 = more output)
##'        indicating whether a warning is given if the required precision
##'        'abstol' has not been reached.
##' @param log logical indicating whether. the logarithmic density is to be computed
##' @param ... see ?pnvmix
##' @author Erik Hintz and Marius Hofert
##' @return numeric vector with the computed probabilities and attributes "error"
##'         (error estimate of the RQMC estimator) and "numiter" (number of iterations)
dnvmixcopula <- function(u, qmix, scale = diag(d), factor = NULL, 
                         control = list(), verbose = FALSE, log = FALSE, ...)
{
   ## Checks 
   if(!is.matrix(u)) u <- rbind(u)
   d <- ncol(u) 
   n <- nrow(u)
   stopifnot(all(u <= 1), all(u >= 0)) 
   ## Result object
   res <- rep(-Inf, n)
   notNA <- rowSums(is.na(u)) == 0 
   not01 <- rowSums( u <= 0 | u >= 1 ) == 0 # rows where no component is <= 0 or >=1
   ## Fill in NAs where needed
   res[!notNA] <- NA
   ## Density is zero outside (0,1)^d 
   res[!not01 & notNA] <- if(log) -Inf else 0 
   u <- u[notNA & not01,, drop = FALSE] # non-missing data inside (0,1)^d
   ## Change accuracy for logdensity in qnvmix() to the one from dnvmix() here
   ## if the user did not provide a different one
   ## (The default for abstol.newton.logdensity is chosen somewhat large for
   ## efficiency reasons as the logdensity there is only needed for Newton)
   names.control <- names(control)
   if(!any(names.control == "newton.logdens.abstol")) {
      ## 'newton.logdens.abstol' was *not* provided:
      control <- get_set_param(control)
      control$newton.logdens.abstol <- control$dnvmix.abstol
   }
   ## Obtain quantiles. Note that qnvmix() takes in and returns a vector
   qu <- qnvmix(as.vector(u), qmix = qmix, control = control,
                verbose = verbose, q.only =  FALSE, ...) # length n*d
   ## log f_{X, scale} (F_{X1}^{-1}(u_{j1}),...,F_X1 ^{-1}(u_{jd})), j = 1,...,n
   num <- dnvmix(matrix(qu$q, ncol = d), qmix = qmix, scale = scale, factor = factor,
                 control = control, verbose = verbose, log = TRUE, ...)# length n
   ## sum_{i=1}^d log f_{X1}( F_{X1}^{-1}(u_{ji})), j = 1,...,n
   ## Note that the log-density values were already calculated by qnvmix()
   denom <- rowSums(matrix(qu$log.density, ncol = d)) # length n
   ## Store results and return
   res[notNA & not01] <- if(!log) exp(num - denom) else num - denom
   res
}


##' @title Distribution Function of a Multivariate Normal Variance Mixture Copula
##' @param upper (n,d) matrix of upper evaluation points. Have to be in (0,1)
##' @param lower (n,d) matrix of lower evaluation points. Have to be in (0,1)
##' @param qmix see ?pnvmix
##' @param scale (d, d)-covariance matrix (scale matrix).
##' @param control see ?get_set_param()
##' @param verbose logical (or integer: 0 = FALSE, 1 = TRUE, 2 = more output)
##'        indicating whether a warning is given if the required precision
##'        'abstol' has not been reached.
##' @author Erik Hintz and Marius Hofert
##' @return numeric vector with the computed probabilities and attributes "error"
##'         (error estimate of the RQMC estimator) and "numiter" (number of iterations)
pnvmixcopula <- function(upper, lower = matrix(0, nrow = n, ncol = d), qmix, 
                      scale = diag(d), control = list(), verbose = FALSE, ...)
{
    ## Most arguments are checked by qnvmix() and pnvmix()
    if(!is.matrix(upper)) upper <- rbind(upper)
    d <- ncol(upper) 
    n <- nrow(upper)
    upper <- pmax( pmin(upper, 1), 0) 
    lower <- pmax( pmin(lower, 1), 0) 
    ## Obtain quantiles. Note that qnvmix() returns a vector
    upper_ <- matrix(qnvmix(as.vector(upper), qmix = qmix, control = control,
                        verbose = verbose, q.only = TRUE, ...), ncol = d)
    lower_ <- if(all(lower == 0)){
       matrix(-Inf, nrow = n, ncol = d) # avoid estimation of the quantile 
    } else {
       matrix(qnvmix(as.vector(lower), qmix = qmix, control = control,
                     verbose = verbose, q.only = TRUE, ...), ncol = d)
    }
    ## Call pnvmix() which handles NA correctly 
    pnvmix(upper_, lower = lower_, qmix = qmix, scale = scale, control = control,
           verbose = verbose, ...)
}


##' @title Random Number Generator for Multivariate Normal Variance Mixtures
##' @param n sample size
##' @param qmix see ?pnvmix
##' @param scale (d, d)-covariance matrix (scale != covariance matrix here)
##' @param factor (d, k)-matrix such that factor %*% t(factor) = scale;
##'        internally determined via chol() (and then an upper triangular
##'        matrix) if not provided
##' @param method character string indicating the method to be used:
##'         - "PRNG":    pure Monte Carlo
##'         - "sobol":   Sobol sequence
##'         - "ghalton": generalized Halton sequence
##'         Note: For the methods "sobol" and "ghalotn", qmix() must be provided
##'         and rmix() is ignored. For the method "PRNG", either qmix() or rmix()
##'         needs to be provided. If both are provided, qmix() is ignored and
##'         rmix() is used.
##' @param skip numeric integer. How many points should be skipped when method='sobol'?
##' @param control see ?get_set_param()
##' @param verbose indicating whether a warning is given if the required precision
##'        in the underlying pnvmix() is not reached.
##' @author Erik Hintz and Marius Hofert
##' @return (n, d)-matrix with NVM(0, scale)-copula samples
rnvmixcopula <- function(n, qmix, scale = diag(2), factor = NULL,
                      method = c("PRNG", "sobol", "ghalton"), skip = 0,
                      control = list(), verbose = FALSE, ...)
{
    d <- dim(scale)[1]
    stopifnot(dim(scale) == c(d,d))
    scale <- cov2cor(scale) # only need correlation matrix

    ## Sample from the nvmix distribution
    sample.nvmix <- rnvmix(n = n, qmix = qmix, scale = scale, factor = factor,
                           method = method, skip = skip, ...)
    ## Apply univariate margins
    ## Need (n,1) matrix as input so that pnvmix() gets the dimension right:
    sample.nvmixcopula <- pnvmix(upper = matrix(sample.nvmix, ncol = 1), qmix = qmix,
                              scale = matrix(1), standardized = TRUE,
                              control = control, verbose = verbose, ...)
    ## Get dimensions correct and return
    matrix(sample.nvmixcopula, ncol = d)
}