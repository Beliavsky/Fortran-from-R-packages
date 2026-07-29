### d/p/rNorm() ################################################################

##' @title Density of the Multivariate Normal Distribution
##' @param x (n, d)-matrix of evaluation points
##' @param loc d-vector (location = mean vector here)
##' @param scale (d, d)-covariance matrix, positive definite (scale = covariance
##'        matrix here)
##' @param factor *lower triangular* factor R of the covariance matrix 'scale'
##'        such that R^T R = 'scale' here (otherwise det(scale) not computed
##'        correctly!)
##' @param log logical indicating whether the logarithmic density is computed
##' @param verbose logical indicating whether a warning is given if the required
##'        precision 'abstol' (see dnvmix()) has not been reached.
##' @param ... additional arguments passed to the underlying dnvmix()
##' @return n-vector of N(loc, scale) density values
##' @author Erik Hintz and Marius Hofert
dNorm <- function(x, loc = rep(0, d), scale = diag(d),
                  factor = NULL, # needs to be triangular!
                  log = FALSE, verbose = TRUE, ...)
{
    if(!is.matrix(x)) x <- rbind(x)
    d <- ncol(x) # for 'loc', 'scale'
    dnvmix(x, qmix = "constant", loc = loc, scale = scale,
           factor = factor, log = log, verbose = verbose, ...)
}


##' @title Distribution Function of the Multivariate Normal Distribution
##' @param upper d-vector of upper evaluation limits
##' @param lower d-vector of lower evaluation limits
##' @param loc d-vector (location = mean vector here)
##' @param scale (d, d)-covariance matrix (scale = covariance matrix here)
##' @param standardized logical indicating whether 'scale' is assumed to be a
##'        correlation matrix; if FALSE (default), 'upper', 'lower' and 'scale'
##'        will be normalized.
##' @param method character string indicating the method to be used:
##'         - "sobol":   Sobol sequence
##'         - "ghalton": generalized Halton sequence
##'         - "prng":    pure Monte Carlo
##' @param precond logical; if TRUE (recommended), variable reordering
##'        similar to Genz and Bretz (2002, pp. 955--956) is performed.
##'        Variable reordering can lead to a significant variance reduction
##'        and decrease in computational time.
##' @param abstol numeric >= 0 providing the absolute precision required.
##'        If abstol = 0, algorithm will run until total number of function
##'        evaluations exceeds fun.eval[2].
##' @param CI.factor Monte Carlo confidence interval multiplier. Algorithm runs
##'        CI.factor * (estimated standard error) < abstol. If CI.factor = 3.3
##'        (default), one can expect the actual absolute error to be less than
##'        abstol in 99.9% of the cases
##' @param fun.eval 2-vector giving the initial function evaluations (in the
##'        first loop; typically powers of 2) and the maximal number of
##'        function evaluations
##' @param B number of randomizations to get error estimates.
##' @param verbose logical indicating whether a warning is given if the required
##'        precision 'abstol' (see dnvmix()) has not been reached.
##' @return numeric vector with the computed probabilities and attributes "error"
##'         (error estimate of the RQMC estimator) and "numiter"
##'         (number of iterations)
##' @author Erik Hintz and Marius Hofert
pNorm <- function(upper, lower = matrix(-Inf, nrow = n, ncol = d),
                  loc = rep(0, d), scale = diag(d), standardized = FALSE,
                  control = list(),
                  verbose = TRUE)
{
    ## Checks (needed to get the default for 'lower' correctly)
    if(!is.matrix(upper)) upper <- rbind(upper) # 1-row matrix if upper is a vector
    n <- nrow(upper) # number of evaluation points
    d <- ncol(upper) # dimension
    pnvmix(upper, lower = lower, qmix = "constant", loc = loc, scale = scale,
           standardized = standardized, control = control, verbose = verbose)
}


##' @title Random Number Generator for the Multivariate Normal Distribution
##' @param n sample size
##' @param loc d-vector (location = mean vector here)
##' @param scale (d, d)-covariance matrix (scale = covariance matrix here)
##' @param factor factor R of the covariance matrix 'scale' with d rows
##'        such that R R^T = 'scale'.
##' @return (n, d)-matrix with N(loc, scale) samples
##' @author Erik Hintz and Marius Hofert
rNorm <- function(n, loc = rep(0, d), scale = diag(2), factor = NULL, # needs to be triangular!
                  method = c("PRNG", "sobol", "ghalton"), skip = 0)
{
   d <- if(!is.null(factor)) { # for 'loc', 'scale'
      nrow(factor <- as.matrix(factor))
   } else {
      nrow(scale <- as.matrix(scale))
   }
   method <- match.arg(method) 
   if(method == "PRNG"){
      ## Provide 'rmix' and no 'qmix' => typically faster
      rnvmix(n, rmix = "constant", loc = loc, scale = scale, factor = factor, 
             method = method, skip = skip)
   } else {
      ## Provide 'qmix' for inversion based methods (needed internally 
      ## even though mixing rv is constant)
      rnvmix(n, qmix = "constant", loc = loc, scale = scale, factor = factor, 
             method = method, skip = skip)
   }
}



##' @title Random Number Generator Multivariate Standard Normal Distribution 
##'        under a Weighted Sum Constraint
##' @param n sample size
##' @param weights d-vector giving weights; all weights must be non-zero
##' @param s either numeric or n-vector of "target sums"      
##' @return (n, d)-matrix of realizations of a standard multivariate normal 
##'         Z=(Z1,..,Zd) conditional on w^T Z = s if 's' is a single number. If
##'         's' is a vector, the i'th row of the return matrix uses the contraint
##'         w^T Z = s_i.
##' @author Erik Hintz 
##' @note Implementation based on Algorithm 1 in Vrins (2018), "Sampling the 
##'       Multivariate Standard Normal Distribution under a Weighted Sum Constraint".

rNorm_sumconstr <- function(n, weights, s, 
                            method = c("PRNG", "sobol", "ghalton"), skip = 0)
{
   ## Basic checks
   stopifnot(n >= 1, all(weights != 0))
   method <- match.arg(method) 
   if(!is.vector(weights)) weights <- as.vector(weights)
   d <- length(weights) # dimension
   method <- match.arg(method)
   ## Turn 's' into a n-vector 
   if(length(s) == 1) s <- rep(s, n) else if(length(s) != n) 
      stop("'s' must be either a single number or a vector of length n")
   norm.w.sq <- sum(weights^2)
   weights.dm1sq  <- weights[1:(d-1)]^2 
   ## Initialize result matrix 'Z.n' 
   Z.n <- matrix(NA, ncol = d, nrow = n)
   ## Compute scale matrix for (Z_1,...,Z_{d-1}) | w^T Z = s (Eq. 3 in Vrins(2018))
   condSig <- -outer(weights.dm1sq, weights.dm1sq)
   diag(condSig) <- weights.dm1sq * (norm.w.sq - weights.dm1sq)
   condSig <- condSig / norm.w.sq
   ## Sample Z_1,...,Z_{d-1} based on condSig
   Z.n[, 1:(d-1)] <- rNorm(n, scale = condSig, method = method, skip = skip) + 
      matrix(weights.dm1sq, ncol = d - 1, nrow = n, byrow = TRUE) * s
   ## Set d'th element correctly to meet constraint
   Z.n[, d] <- s - rowSums(Z.n[, 1:(d-1), drop = FALSE])
   Z.n <- Z.n * matrix(1/weights, ncol = d, nrow = n, byrow = TRUE) 
   ## Return
   Z.n 
}




##' @title Fitting the Parameters of a Multivariate Normal Distribution
##' @param x (n,d) data matrix
##' @return list with components 'loc' (sample mean) and 'scale'
##'         (sample covariance matrix)
##' @author Marius Hofert
fitNorm <- function(x)
{
    if(!is.matrix(x))
        x <- rbind(x)
    list(loc = colMeans(x), scale = cov(x))
}
