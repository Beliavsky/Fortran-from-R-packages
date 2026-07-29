### GPD-based tail distribution (POT) ##########################################

##' @title Density of the GPD-Based Tail Distribution (POT)
##' @param x evaluation points
##' @param threshold threshold u
##' @param p.exceed probability of exceeding the threshold u;
##'        for Smith estimator, this would be mean(x > threshold) for
##'        'x' being the data.
##' @param shape parameter xi
##' @param scale parameter beta
##' @param log logical indicating whether the log density is computed
##' @return f(x) = bar{F}_n(u) * f_GPD(x - threshold)
##' @author Marius Hofert
dGPDtail <- function(x, threshold, p.exceed, shape, scale, log = FALSE)
{
    stopifnot(0 <= p.exceed, p.exceed <= 1, scale > 0) # don't need to check x > threshold as this is done by dGPD()
    res <- log(p.exceed) + dGPD(x - threshold, shape = shape, scale = scale, log = TRUE)
    if(log) res else exp(res)
}

##' @title GPD-Based Tail Distribution (POT Method)
##' @param q quantile (vectorized)
##' @param threshold threshold u
##' @param p.exceed probability of exceeding the threshold u;
##'        for Smith estimator, this would be mean(x > threshold) for
##'        'x' being the data.
##' @param shape parameter xi
##' @param scale parameter beta
##' @param lower.tail logical indicating whether lower/upper tail is used
##' @param log.p logical indicating whether probabilities are given as log()
##' @return F(q) = 1-bar{F}_n(u) (1-F_GPD(q-u))
##' @author Marius Hofert
pGPDtail <- function(q, threshold, p.exceed, shape, scale, lower.tail = TRUE, log.p = FALSE)
{
    stopifnot(q >= threshold, 0 <= p.exceed, p.exceed <= 1, scale > 0)
    res <- p.exceed * pGPD(q-threshold, shape = shape, scale = scale, lower.tail = FALSE) # = bar(F(q))
    if(lower.tail) {
        if(log.p) {
            log1p(-res) # = log(F(q))
        } else {
            1-res # F(q)
        }
    } else {
        if(log.p) {
            log(res) # = log(bar(F(q))); no better log available (it seems)
        } else {
            res
        }
    }
}

##' @title Quantile function of the GPD-Based Tail Distribution (POT Method)
##' @param p probability (vectorized)
##' @param threshold threshold u
##' @param p.exceed probability of exceeding the threshold u;
##'        for Smith estimator, this would be mean(x > threshold) for
##'        'x' being the data.
##' @param shape parameter xi
##' @param scale parameter beta
##' @param lower.tail logical indicating whether lower/upper tail is used
##' @param log.p logical indicating whether probabilities are given as log()
##' @return F^{-1}(q) = u + F_GPD^{-1}(1-(1-p)/bar{F}_n(u))
##' @author Marius Hofert
qGPDtail <- function(p, threshold, p.exceed, shape, scale, lower.tail = TRUE, log.p = FALSE)
{
    stopifnot(0 <= p.exceed, p.exceed <= 1, scale > 0)
    p <- if(log.p) pmin(p, 0) else pmin(pmax(p, 0), 1) # correctly extend
    if(lower.tail) {
        if(log.p) {
            threshold + qGPD(1-expm1(-p)/p.exceed, shape = shape, scale = scale)
        } else {
            threshold + qGPD(1-(1-p)/p.exceed, shape = shape, scale = scale)
        }
    } else {
        if(log.p) {
            threshold + qGPD(1-exp(p)/p.exceed, shape = shape, scale = scale)
        } else {
            threshold + qGPD(1-p/p.exceed, shape = shape, scale = scale)
        }
    }
}

##' @title Generating random variates from a GPD-Based Tail Distribution (POT Method)
##' @param n sample size n
##' @param threshold threshold u
##' @param p.exceed probability of exceeding the threshold u;
##'        for Smith estimator, this would be mean(x > threshold) for
##'        'x' being the data.
##' @param shape parameter xi
##' @param scale parameter beta
##' @return n-vector containing random variates
##' @author Marius Hofert
rGPDtail <- function(n, threshold, p.exceed, shape, scale)
    qGPDtail(runif(n), threshold = threshold, p.exceed = p.exceed,
             shape = shape, scale = scale)
