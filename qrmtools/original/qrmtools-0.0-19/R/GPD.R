### GPD(shape, scale) distribution #############################################

##' @title Density of the GPD(shape, scale) distribution
##' @param x evaluation points
##' @param shape parameter xi
##' @param scale parameter beta
##' @param log logical indicating whether the log density is computed
##' @return density of the GPD(shape, scale) distribution
##' @author Marius Hofert
dGPD <- function(x, shape, scale, log = FALSE)
{
    l <- length(x)
    if(scale <= 0)
        return(rep(if(log) -Inf else 0, l)) # for logLik_GPD()
    if(shape == 0) { # shape == 0
        res <- -x/scale-log(scale)
    } else { # shape != 0
        ## Note: If shape < 0, the support is [0, -scale/shape]
        res <- rep(-Inf, l) # correctly extend log-density
        ii <- if(shape > 0) 0 <= x else 0 <= x & x < -scale/shape # those indices for which density is positive
        res[ii] <- -(1/shape + 1) * log1p(shape * x[ii] / scale) - log(scale)
    }
    if(log) res else exp(res)
}

##' @title Distribution function of the GPD(shape, scale) distribution
##' @param q quantile (vectorized)
##' @param shape parameter xi
##' @param scale parameter beta
##' @param lower.tail logical indicating whether lower/upper tail is used
##' @param log.p logical indicating whether probabilities are given as log()
##' @return distribution function of the GPD(shape, scale) distribution
##' @author Marius Hofert
pGPD <- function(q, shape, scale, lower.tail = TRUE, log.p = FALSE)
{
    stopifnot(scale > 0)
    ## Note: The following extension doesn't work as there can be cases (2022-02-26)
    ##       when q = -scale/shape but 1 + q * shape/scale < 0 at which
    ##       point the below code for shape != 0 fails.
    ## q <- if(shape >= 0) pmax(q, 0) else pmin(pmax(q, 0), -scale/shape)
    res <- numeric(length(q))
    if(shape >= 0) {
        ii <- q <= 0 # include boundary case here to not run into numerical problems
        res[ii] <- 0
    } else {
        i <- q <= 0 # include boundary case here to not run into numerical problems
        j <- q >= -scale/shape # include boundary case here to not run into numerical problems
        res[i] <- 0
        res[j] <- 1
        ii <- i | j
    }
    res[!ii] <- if(shape == 0) { # shape == 0
                    if(lower.tail) {
                        if(log.p) {
                            log1p(-exp(-q[!ii]/scale))
                        } else {
                            1-exp(-q[!ii]/scale)
                        }
                    } else {
                        if(log.p) {
                            -q[!ii]/scale
                        } else {
                            exp(-q[!ii]/scale)
                        }
                    }
                } else { # shape != 0
                    if(lower.tail) {
                        if(log.p) {
                            log1p(-(1+shape*q[!ii]/scale)^(-1/shape))
                        } else {
                            1-(1+shape*q[!ii]/scale)^(-1/shape)
                        }
                    } else {
                        if(log.p) {
                            -log1p(shape*q[!ii]/scale)/shape
                        } else {
                            (1+shape*q[!ii]/scale)^(-1/shape)
                        }
                    }
                }
    res
}

##' @title Quantile function of GPD(shape, scale)
##' @param p probability (vectorized)
##' @param shape parameter xi
##' @param scale parameter beta
##' @param lower.tail logical indicating whether lower/upper tail is used
##' @param log.p logical indicating whether probabilities are given as log()
##' @return quantile function of the GPD(shape, scale) distribution
##' @author Marius Hofert
qGPD <- function(p, shape, scale, lower.tail = TRUE, log.p = FALSE)
{
    stopifnot(scale > 0)
    p <- if(log.p) pmin(p, 0) else pmin(pmax(p, 0), 1) # correctly extend
    if(shape == 0) { # shape == 0
        if(lower.tail)
            if(log.p) (-scale)*log1p(-exp(p)) else (-scale)*log1p(-p)
        else if(log.p) (-scale)*p else (-scale)*log(p)
    } else { # shape != 0
        if(lower.tail)
            if(log.p) (scale/shape)*((-expm1(p))^(-shape)-1) else (scale/shape)*((1-p)^(-shape)-1)
        else if(log.p) (scale/shape)*expm1(-shape*p) else (scale/shape)*(p^(-shape)-1)
    }
}

##' @title Generating random variates from a GPD(shape, scale) distribution
##' @param n sample size n
##' @param shape parameter xi
##' @param scale parameter beta
##' @return n-vector containing GPD(shape, scale) random variates
##' @author Marius Hofert
rGPD <- function(n, shape, scale)
    qGPD(runif(n), shape = shape, scale = scale)


### Par(shape, scale) = Par(theta, kappa) = GPD(1/theta, kappa/theta), theta > 0 distribution

## Note: - Hard-coded here to be vectorized in the main argument and theta
##       - F(x) = 1 - (1+x/kappa)^{-theta}, theta > 0, kappa > 0, x >= 0
##       - E[X] = kappa / (theta-1) for all theta > 1 (see McNeil, Frey, Embrechts (2015))
##       - Var[X] = theta * kappa^2 / ((theta-2)(theta-1)^2) for all theta > 2 (see McNeil, Frey, Embrechts (2015))

##' @title Density of the Par(shape, scale) distribution
##' @param x evaluation points
##' @param shape parameter theta
##' @param scale parameter kappa
##' @param log logical indicating whether the log density is computed
##' @return density of the Par(shape, scale) distribution
##' @author Marius Hofert
dPar <- function(x, shape, scale = 1, log = FALSE)
{
    stopifnot(shape > 0, scale > 0)
    if(log) log(shape/scale) + (shape+1) * log(scale/(scale+x)) else (shape/scale)*(scale/(scale+x))^(shape+1)
}

##' @title Distribution function of the Par(shape, scale) distribution
##' @param q quantile
##' @param shape parameter theta
##' @param scale parameter kappa
##' @param lower.tail logical indicating whether lower/upper tail is used
##' @param log.p logical indicating whether probabilities are given as log()
##' @return distribution function of the Par(shape, scale) distribution
##' @author Marius Hofert
pPar <- function(q, shape, scale = 1, lower.tail = TRUE, log.p = FALSE)
{
    stopifnot(shape > 0, scale > 0)
    if(lower.tail) {
        if(log.p) log(1-(scale/(scale+q))^shape) else 1-(scale/(scale+q))^shape
    } else if(log.p) shape*(log(scale)-log(scale+q)) else (scale/(scale+q))^shape
}

##' @title Quantile function of the Par(shape, scale) distribution
##' @param p probability
##' @param shape parameter theta
##' @param scale parameter kappa
##' @param lower.tail logical indicating whether lower/upper tail is used
##' @param log.p logical indicating whether probabilities are given as log()
##' @return quantile function of the Par(shape, scale) distribution
##' @author Marius Hofert
qPar <- function(p, shape, scale = 1, lower.tail = TRUE, log.p = FALSE)
{
    stopifnot(0 <= p, p <= 1, shape > 0, scale > 0)
    if(lower.tail) {
        if(log.p) scale * ((-expm1(p))^(-1/shape)-1) else scale * ((1-p)^(-1/shape)-1)
    } else if(log.p) scale * expm1(-p/shape) else scale * (p^(-1/shape)-1)
}

##' @title Generating random variates from a Pareto(shape, scale) distribution
##' @param n sample size n
##' @param shape parameter theta
##' @param scale parameter kappa
##' @return n-vector containing Pareto(shape, scale) random variates
##' @author Marius Hofert
rPar <- function(n, shape, scale = 1)
{
    stopifnot(shape > 0, scale > 0)
    qPar(runif(n), shape = shape, scale = scale)
}

##' @title Primitive of the Par(shape, scale) survival function
##' @param q quantile
##' @param shape parameter theta
##' @param scale parameter kappa
##' @return \int\bar{F}(x) dx
##' @author Marius Hofert
bar_pPar_primitive <- function(q, shape, scale = 1)
{
    stopifnot(shape > 0, scale > 0)
    if(shape == 1) scale*log(scale+q) else (scale/(1-shape)) * (scale/(scale+q))^(shape-1)
}

