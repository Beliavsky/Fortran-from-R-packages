### Allocation #################################################################

##' @title Euler Allocations for Elliptical Distributions
##' @param total total to be allocated (typically the risk measure of the
##'        sum of the underlying loss random variables)
##' @param loc location vector of the elliptical distribution of the loss
##'        random vector
##' @param Sigma scale (covariance) matrix of the elliptical distribution of
##'        the loss random vector
##' @return allocation according to the Euler principle
##' @author Marius Hofert and Takaaki Koike
##' @note 1) See MFE (2015, Corollary 8.43) for bm(L) ~ E_d(bm{0}, Sigma, psi)
##'          and positive-homogeneous and law-invariant risk measures.
##'          After summing over all i, you get AC (total) / AC_j = sum(Sigma) /
##'          rowSums(Sigma), or, equivalently, rowSums(Sigma) / sum(Sigma)
##'          = AC_j / AC. Since AC = total, we obtain AC_j = (rowSums(Sigma) /
##'          sum(Sigma)) * total.
##'       2) If, additionally, the risk measure is also translation invariant,
##'          MFE (2015, Theorem 8.28 (1)) implies that r_varrho(bm(lambda))
##'          gets an additional summand bm(lambda)^T bm(mu). When differentiated
##'          w.r.t. lambda_i, this implies a summand of mu_i, so tilde(AC)_j
##'          := AC_j - mu_j satisfies 1) and thus tilde(AC)_j = (rowSums(Sigma) /
##'          sum(Sigma)) * total and thus AC_j = mu_j + (rowSums(Sigma) /
##'          sum(Sigma)) * total
alloc_ellip <- function(total, loc, scale)
{
    d <- ncol(scale)
    if(length(loc) == 1) loc <- rep(loc, d)
    stopifnot(nrow(scale) == d, length(loc) == d)
    rs <- rowSums(scale)
    loc + (rs / sum(rs)) * total # Euler allocation
}

##' @title Sub-sampling based on a Risk Measure of the Sum
##' @param x (n, d)-data matrix
##' @param level confidence level(s)
##' @param risk.measure character string or function specifying the risk measure
##'        computed from the row sums of x based on the given level(s) in order
##'        to determine the conditioning region.
##' @param ... additional arguments passed to 'risk.measure'
##' @return sub-sample of x that satisfies that its row sums are in the
##'         conditioning region specified by 'risk.measure' computed from the row
##'         sums of x and the given 'level'
##' @author Marius Hofert
##' @note This function could at some point get an argument 'type' for using,
##'       say, the row maximum instead of the row sum to determine the
##'       conditioning region.
conditioning <- function(x, level, risk.measure = "VaR_np", ...)
{
    ## Basics
    if(!is.matrix(x))
        x <- as.matrix(x)
    if(length(level) == 1) level <- c(level, 1)
    stopifnot(length(level) == 2, 0 <= level, level <= 1)
    is.function.rm <- is.function(risk.measure)
    if(!is.function.rm) {
        stopifnot(is.character(risk.measure), existsFunction(risk.measure)) # check
        ## Create an expression of the function call and evaluate that below
        expr <- as.call(c(as.name(risk.measure), # 'unquote' string
                          quote(x), # includes 'x' as first argument (which exists in here)
                          quote(level.), # includes placeholder to be substituted below
                          as.expression(list(...)))) # includes '...' (which exists in here)
    }

    ## Estimate risk.measure(S) for the two confidence levels
    ## Note: We could also call risk.measure(S, level = level, ...) directly
    ##       but that assumes the provided risk.measure() is vectorized in 'level'
    S <- rowSums(x) # row sums
    rm.level.S <- if(is.function.rm) {
                      sapply(level, function(l) risk.measure(S, level = l, ...))
                  } else {
                      sapply(level, function(l) eval(expr, list(level. = l)))
                  }
    if(length(rm.level.S) != 2) # sanity check
        stop("The output of the evaluated 'risk.measure' does not have length 2.")

    ## Return sub-sample
    x[(rm.level.S[1] < S) & (S <= rm.level.S[2]),]
}

##' @title Nonparameteric Allocation
##' @param x see ?conditioning
##' @param level see ?conditioning
##' @param risk.measure see ?conditioning
##' @param include.conditional logical indicating whether the conditional data set
##'        is returned, too
##' @param ... see ?conditioning
##' @return list with components:
##'         "allocation": nonparametric allocation estimate
##'         "SE": standard error of the allocation amounts
##'         "n": sample size of the conditional sample based on which the allocation
##'              is estimated
##'         "conditional": the underlying conditional sample of X given that S lies in the
##'                        region defined by the provided risk measure; omitted if
##'                        include.conditional = FALSE
##' @author Marius Hofert
##' @note For an iid sample x from X (assumed here), note that the SE is the standard
##'       deviation of the sample mean (since Var(<sample mean>) = Var(X_1)/n
##'       => sd(<sample mean>) = sigma/sqrt{n} = SE), so a measure of how well the
##'       estimated allocation approximates the true one.
alloc_np <- function(x, level, risk.measure = "VaR_np", include.conditional = FALSE, ...)
{
    z <- conditioning(x, level = level, risk.measure = risk.measure, ...) # Z = X | varrho(S) in region
    if(!is.matrix(z)) z <- as.matrix(z)
    res <- list(allocation = colMeans(z), # estimated allocation
                SE = apply(z, 2, sd) / sqrt(nrow(z)), # estimated standard errors per allocated amount
                n = nrow(z))
    if(include.conditional) c(res, conditional = z) else res
}
