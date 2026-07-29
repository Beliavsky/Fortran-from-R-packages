### GEV(shape, loc, scale) distribution ########################################

##' @title Density of the GEV(shape, loc, scale) distribution
##' @param x evaluation points
##' @param shape parameter xi (real)
##' @param loc parameter mu (real)
##' @param scale parameter sigma (also real here; density is 0 if sigma <= 0)
##' @param log logical indicating whether the log density is computed
##' @return density of the GEV(xi, mu, sigma) distribution
##' @author Marius Hofert
dGEV <- function(x, shape, loc = 0, scale = 1, log = FALSE)
{
    l <- length(x)
    if(scale <= 0)
        return(rep(if(log) -Inf else 0, l)) # for logLik_GEV()
    y <- (x - loc) / scale # scale > 0
    if(shape == 0) { # shape == 0
        res <- -log(scale) - (y + exp(-y))
        ## Note: - for loc >> x, y << 0 => exp(-y) = Inf => density correctly 0 then
        ##       - y should be >= -log(.Machine$double.xmax) for exp(-y) < Inf
    } else { # shape != 0
        res <- rep(-Inf, l) # correctly extend log-density
        sy <- shape * y
        ii <- 1 + sy > 0 # those indices for which density is positive
        res[ii] <- -log(scale) + (-1/shape - 1) * log1p(sy[ii]) - (1 + sy[ii])^(-1/shape)
    }
    if(log) res else exp(res)
}

##' @title Distribution function of the GEV(shape, loc, scale) distribution (vectorized in q)
##' @param q quantile
##' @param shape parameter xi
##' @param loc parameter mu
##' @param scale parameter sigma
##' @param lower.tail logical indicating whether lower/upper tail is used
##' @param log.p logical indicating whether probabilities are given as log()
##' @return distribution function of the GEV(shape, loc, scale) distribution
##' @author Marius Hofert
pGEV <- function(q, shape, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
{
    stopifnot(scale > 0)
    y <- (q - loc) / scale
    if(shape == 0) { # shape == 0
        if(lower.tail)
            if(log.p) -exp(-y) else exp(-exp(-y))
        else if(log.p) log1p(-exp(-exp(-y))) else 1-exp(-exp(-y))
    } else { # shape != 0
        sy <- pmax(shape*y, -1) # see dGEV()
        if(lower.tail) {
            if(log.p) {
                -(1+sy)^(-1/shape) # log H
            } else {
                exp(-(1+sy)^(-1/shape)) # H
            }
        } else {
            if(log.p) {
                log1p(-exp(-(1+sy)^(-1/shape))) # log(bar{H})
            } else {
                1-exp(-(1+sy)^(-1/shape)) # bar{H}
            }
        }
    }
}

##' @title Quantile function of the GEV(shape, loc, scale) distribution (vectorized in p)
##' @param p probability
##' @param shape parameter xi
##' @param loc parameter mu
##' @param scale parameter sigma
##' @param lower.tail logical indicating whether lower/upper tail is used
##' @param log.p logical indicating whether probabilities are given as log()
##' @return quantile function of the GEV(shape, loc, scale) distribution
##' @author Marius Hofert
qGEV <- function(p, shape, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
{
    stopifnot(scale > 0)
    p <- if(log.p) pmin(p, 0) else pmin(pmax(p, 0), 1) # correctly extend
    if(shape == 0) { # shape == 0
        res <- if(lower.tail)
                   if(log.p) -log(-p) else -log(-log(p))
               else if(log.p) -log(-log1p(-exp(p))) else -log(-log1p(-p))
    } else { # shape != 0
          res <- if(lower.tail)
                     if(log.p) ((-p)^(-shape)-1)/shape else ((-log(p))^(-shape)-1)/shape
                 else if(log.p) ((-log1p(-exp(p)))^(-shape)-1)/shape else ((-log1p(-p))^(-shape)-1)/shape
          res <- if(shape < 0) pmin(res, -1/shape) else pmax(res, -1/shape)
    }
    loc + scale * res
}

##' @title Generating random numbers from a GEV(shape, loc, scale) distribution
##' @param n sample size n
##' @param shape parameter xi
##' @param loc parameter mu
##' @param scale parameter sigma
##' @return n-vector containing GEV(shape, loc, scale) random variates
##' @author Marius Hofert
rGEV <- function(n, shape, loc = 0, scale = 1)
    qGEV(runif(n), shape = shape, loc = loc, scale = scale)


