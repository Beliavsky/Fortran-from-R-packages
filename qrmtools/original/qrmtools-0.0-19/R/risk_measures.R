### Computing risk measures ####################################################

### 1 Value-at-Risk ############################################################

##' @title Nonparametric VaR Estimator
##' @param x vector of losses
##' @param level confidence level alpha
##' @param names see ?quantile
##' @param type 'type' used (1 = inverse of empirical df); see ?quantile
##' @param ... see ?quantile
##' @return nonparametric VaR_alpha estimate
##' @author Marius Hofert
##' @note We use the conservative type = 1 here as for sufficiently large alpha,
##'       type = 7 (quantile()'s default) would interpolate between the two
##'       largest losses (to be continuous) and thus return a(n even) smaller
##'       VaR_alpha estimate.
VaR_np <- function(x, level, names = FALSE, type = 1, ...)
    quantile(if(is.matrix(x)) rowSums(x) else x, # compute VaR of the sum if a matrix is provided
             probs = level, names = names, type = type, ...) # vectorized in level

##' @title Value-at-Risk for the t Distribution
##' @param level confidence level alpha
##' @param loc location mu
##' @param scale scale sigma
##' @param df degrees of freedom; Inf for the normal distribution
##' @return Value-at-Risk
##' @author Marius Hofert
VaR_t <- function(level, loc = 0, scale = 1, df = Inf)
{
    stopifnot(0 <= level, level <= 1, scale > 0, df > 0)
    loc + scale * if(identical(df, Inf)) qnorm(level) else qt(level, df = df)
}

##' @title Value-at-Risk for the Standardized t Distribution
##' @param level confidence level alpha
##' @param df degrees of freedom; Inf for the normal distribution
##' @return Value-at-Risk
##' @author Marius Hofert
VaR_t01 <- function(level, df = Inf)
{
    if(df <= 2)
        stop("Standardized t distribution requires df > 2")
    c <- if(identical(df, Inf)) 1 else sqrt((df-2)/df)
    c * VaR_t(level, df = df)
}

##' @title Value-at-Risk for the GPD
##' @param level confidence level alpha
##' @param shape parameter xi
##' @param scale parameter beta
##' @return Value-at-Risk
##' @author Marius Hofert
VaR_GPD <- function(level, shape, scale)
    qGPD(level, shape = shape, scale = scale)

##' @title Value-at-Risk for the Pareto Distribution
##' @param level confidence level alpha
##' @param shape parameter theta
##' @param scale parameter kappa
##' @return Value-at-Risk
##' @author Marius Hofert
VaR_Par <- function(level, shape, scale = 1)
    qPar(level, shape = shape, scale = scale)

##' @title Semi-parametric VaR Estimator in the POT Method
##' @param level confidence level alpha
##' @param threshold threshold u
##' @param p.exceed exceedance probability; typically mean(x > threshold)
##'        for x being the original data
##' @param shape parameter xi
##' @param scale parameter beta
##' @return Value-at-Risk
##' @author Marius Hofert
VaR_GPDtail <- function(level, threshold, p.exceed, shape, scale)
{
    stopifnot(0 <= p.exceed, p.exceed <= level, level <= 1, scale > 0)
    ## equals threshold + (scale/shape) * (((1-level)/p.exceed)^(-shape) - 1)
    qGPDtail(level, threshold = threshold, p.exceed = p.exceed,
             shape = shape, scale = scale)
}


### 2 Expected shortfall #######################################################

##' @title Nonparametric Expected Shortfall Estimator
##' @param x vector of losses
##' @param level confidence level alpha
##' @param method method
##' @param verbose logical indicating whether verbose output is provided in
##'        case the mean is taken over (too) few losses
##' @param ... additional arguments passed to the underlying VaR_np()
##' @return nonparametric ES_{alpha} estimate (derived under the assumption of continuity)
##' @author Marius Hofert
##' @note - Vectorized in level
##'       - ">" : Mathematically correct for continuous and discrete dfs, but
##'               produces NaN for level > (n-1)/n (=> F^-(level) = x_{(n)} but
##'               there is no loss strictly beyond x_{(n)})
##'         ">=": mean() will always include the largest loss (so no NaN appears),
##'               but might be computed just based on this one loss.
ES_np <- function(x, level, method = c(">", ">="), verbose = FALSE, ...)
{
    if(is.matrix(x)) x <- rowSums(x) # compute ES of the sum if a matrix is provided
    stopifnot(0 < level, level < 1)
    VaR <- VaR_np(x, level = level, ...) # length(level)-vector
    method <- match.arg(method)
    vapply(VaR, function(v) { # v = VaR value for one level
        ind <- if(method == ">") x > v else x >= v
        if(verbose) {
            num <- sum(ind)
            if(num == 0) {
                warning("No loss ",method," VaR; NaN returned instead")
            } else if(num == 1){
                warning("Only ",num," loss ",method," VaR")
            } else if(num <= 5) {
                warning("Only ",num," losses ",method," VaR")
            }
        }
        mean(x[ind]) # mean over only those losses >(=) VaR (= E(L; L > VaR(L)) / (1-alpha))
    }, NA_real_)
}

##' @title Expected Shortfall for the t Distribution
##' @param level confidence level alpha
##' @param loc location mu
##' @param scale scale sigma
##' @param df degrees of freedom; Inf for the normal distribution
##' @return Expected shortfall
##' @author Marius Hofert
ES_t <- function(level, loc = 0, scale = 1, df = Inf)
{
    stopifnot(0 <= level, level <= 1, scale > 0, df > 0)
    loc + (scale/(1-level)) * if(identical(df, Inf)) dnorm(qnorm(level)) else
    dt(qt(level, df = df), df = df) * (df + qt(level, df = df)^2) / (df-1)
}

##' @title Expected Shortfall for Standardized t Distribution
##' @param level confidence level alpha
##' @param df degrees of freedom; Inf for the normal distribution
##' @return Expected shortfall
##' @author Marius Hofert
##' @note Because the quantile function scales with c, so does ES
ES_t01 <- function(level, df = Inf)
{
    if(df <= 2)
        stop("Standardized t distribution requires df > 2")
    c <- if(identical(df, Inf)) 1 else sqrt((df-2)/df)
    c * ES_t(level, df = df)
}

##' @title Expected Shortfall for the GPD
##' @param level confidence level alpha
##' @param shape parameter xi
##' @param scale parameter beta
##' @return Expected shortfall
##' @author Marius Hofert
ES_GPD <- function(level, shape, scale)
{
    stopifnot(0 <= level, level <= 1, shape < 1, scale > 0)
    VaR <- VaR_GPD(level, shape = shape, scale = scale)
    VaR + (scale + shape * VaR)/(1-shape) # VaR + mean excess function at VaR
}

##' @title Expected Shortfall for the Pareto Distribution
##' @param level confidence level alpha
##' @param shape parameter theta
##' @param scale parameter kappa
##' @return Expected shortfall
##' @author Marius Hofert
ES_Par <- function(level, shape, scale = 1)
{
    stopifnot(0 <= level, level <= 1, shape > 1, scale > 0)
    scale * ((shape / (shape-1)) * (1-level)^(-1/shape) - 1)
}

##' @title Semi-parametric ES Estimator in the POT Method
##' @param level confidence level alpha
##' @param threshold threshold u
##' @param p.exceed exceedance probability; typically mean(x > threshold)
##'        for x being the original data
##' @param shape parameter xi
##' @param scale parameter beta
##' @return Expected shortfall
##' @author Marius Hofert
ES_GPDtail <- function(level, threshold, p.exceed, shape, scale)
{
    VaR <- VaR_GPDtail(level, threshold = threshold, p.exceed = p.exceed,
                       shape = shape, scale = scale) # does checks
    res <- (VaR + scale - shape * threshold) / (1 - shape)
    res[shape >= 1] <- Inf
    res
}


### 3 Range value-at-risk ######################################################

##' @title Nonparametric Range VaR Estimator
##' @param x vector of losses
##' @param level lower and upper confidence levels; if length(level) == 1, the
##'        upper one is taken as 1
##' @param ... additional arguments passed to the underlying VaR_np()
##' @return nonparametric RVaR_{alpha, beta} estimate
##' @author Marius Hofert
RVaR_np <- function(x, level, ...)
{
    ## Basics
    if(is.matrix(x)) x <- rowSums(x) # compute RVaR of the sum if a matrix is provided
    len <- length(level)
    stopifnot(len == 1 || len == 2, 0 < level, level <= 1, diff(level) > 0)
    if(len == 1) level <- c(level , 1)

    ## Compute nonparametric VaR estimates for both confidence levels
    VaR <- VaR_np(x, level = level, ...)

    ## Estimate RVaR
    mean(x[VaR[1] < x & x <= VaR[2]])
}


### 4 Geometric risk measures ##################################################

##' @title Objective Function E(Lambda_alpha(X-c))
##' @param x running variable ("c")
##' @param X (n, d)-matrix containing random sample from the underlying
##'        multivariate distribution H
##' @param level d-vector of confidence levels
##' @param measure character string specifying the risk measure
##' @return value of E(Lambda_alpha(X-c))
##' @author Marius Hofert
ELambda <- function(x, X, level, measure = c("gExp", "gVaR"))
{
    X. <- sweep(X, MARGIN = 2, STATS = x, FUN = "-")
    norm <- sqrt(rowSums(X.^2))
    measure <- match.arg(measure)
    switch(measure,
    "gExp" = {
        mean(0.5 * norm * (norm + rowSums(X. * rep(level, each = nrow(X.)))))
    },
    "gVaR" = {
        mean(0.5 *        (norm + rowSums(X. * rep(level, each = nrow(X.)))))
    },
    stop("Wrong 'measure'"))
}

##' @title Multivariate Geometric VaR
##' @param x (n, d)-matrix containing random sample from the underlying
##'        multivariate distribution H
##' @param level d-vector of confidence levels; can also be a (N, d)-matrix where
##'        N is the number of 'level' vectors of dimension d each
##' @param start d-vector of initial values for the minimization of ELambda()
##' @param method method for optim()
##' @param ... additional arguments passed to the underlying optim()
##' @return value of optim() containing the geometric VaR in the component 'par'
##' @author Marius Hofert
##' @note if length(level) == 1, 'lower' and 'upper' need to be provided for
##'       method "Brent"
gVaR <- function(x, level, start = colMeans(x),
                 method = if(length(level) == 1) "Brent" else "Nelder-Mead", ...)
{
    if(!is.matrix(x)) x <- rbind(x)
    if(is.matrix(level)) {
        apply(level, 1, function(a) optim(par = start, fn = ELambda, X = x, level = a,
                                          measure = "gVaR", method = method, ...))
    } else { # level is a single vector (for d > 1) or number (for d = 1)
        optim(par = start, fn = ELambda, X = x, level = level, measure = "gVaR",
              method = method, ...)
    }
}

##' @title Multivariate Geometric Expectile
##' @param x (n, d)-matrix containing random sample from the underlying
##'        multivariate distribution H
##' @param level d-vector of confidence levels; can also be a (N, d)-matrix where
##'        N is the number of level vectors of dimension d each
##' @param start d-vector of initial values for the minimization of ELambda()
##' @param method method for optim()
##' @param ... additional arguments passed to the underlying optim()
##' @return value of optim() containing the geometric Exp in the component 'par'
##' @author Marius Hofert
##' @note if length(level) == 1, 'lower' and 'upper' need to be provided for
##'       method "Brent"
gEX <- function(x, level, start = colMeans(x),
                 method = if(length(level) == 1) "Brent" else "Nelder-Mead", ...)
{
    if(!is.matrix(x)) x <- rbind(x)
    if(is.matrix(level)) {
        apply(level, 1, function(a) optim(par = start, fn = ELambda, X = x, level = a,
                                          measure = "gExp", method = method, ...))
    } else { # level is a single vector (for d > 1) or number (for d = 1)
        optim(par = start, fn = ELambda, X = x, level = level, measure = "gExp",
              method = method, ...)
    }
}
