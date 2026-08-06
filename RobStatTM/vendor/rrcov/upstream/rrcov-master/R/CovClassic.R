## Just a shortcut and for competability with older versions
##  same as CovClassic()
##
Cov <- function(x, unbiased = TRUE)
{
    return (CovClassic(x, unbiased))
}

CovClassic <- function(x, unbiased = TRUE)
{
    if(is.data.frame(x))
        x <- data.matrix(x)
    else if(!is.matrix(x))
        x <- matrix(x, length(x), 1,
            dimnames = list(names(x), deparse(substitute(x))))

    ## drop all rows with missing values (!!) :
    na.x <- !is.finite(x %*% rep(1, ncol(x)))
    ok <- !na.x
    x <- x[ok, , drop = FALSE]
    dx <- dim(x)
    if(!length(dx))
        stop("All observations have missing values!")

    dimn <- dimnames(x)
    n <- dx[1]
    p <- dx[2]

##    if(n < p)
##        stop("Need at least p=(number of variables) observations ")

    if(n <= 0)
        stop("All observations have missing values!")

    call = match.call()
    method = "Classical Estimator."
    ans <- cov.wt(x)

    nobs <- nrow(x)
    if(!unbiased)
        ans$cov <- (ans$cov * (nobs-1))/nobs

    new("CovClassic", call=call, cov=ans$cov, center=ans$center, n.obs=ans$n.obs,
        method=method, X=x)
}

## VT::17.06.2008
##setMethod("plot", "CovClassic", function(x, y="missing",
setMethod("plot", signature(x="CovClassic", y="missing"), function(x, y="missing",
                                which=c("all", "distance", "qqchi2", "tolEllipsePlot", "screeplot", "pairs"),
                                ask = (which=="all" && dev.interactive(TRUE)),
                                cutoff,
                                id.n,
                                tol = 1e-7, ...)
{
    data <- getData(x)
    ##  parameters and preconditions
    if(is.vector(data) || is.matrix(data)) {
        if(!is.numeric(data))
            stop( "x is not a numeric dataframe or matrix.")
    } else if(is.data.frame(data)) {
        if(!all(sapply(data,data.class) == "numeric"))
            stop( "x is not a numeric dataframe or matrix.")
    }

    n <- dim(data)[1]
    p <- dim(data)[2]

    if(length(getCenter(x))  == 0 ||  length(getCov(x)) == 0)
        stop( "Invalid object: attributes center and cov missing!")

    if(length(getCenter(x))  != p)
        stop( "Data set and provided center have different dimensions!")

    ## Check for singularity of the cov matrix
    if(isSingular(x))
        stop("The covariance matrix is singular!")

    if(missing(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))

    if(!missing(id.n) && !is.null(id.n)) {
        id.n <- as.integer(id.n)
        if(id.n < 0 || id.n > n)
            stop(sQuote("id.n")," must be in {1,..,",n,"}")
    }

    ccov <- x
    rd <- md <- sqrt(getDistance(x))
    which <- match.arg(which)

    op <- if (ask) par(ask = TRUE) else list()
    on.exit(par(op))

    ## index plot of mahalanobis distances
    if(which == "all" || which == "distance") {
        .mydistplot(md, cutoff, classic=TRUE, id.n=id.n, ...)
    }

    ## qq-plot of the mahalanobis distances versus the
    ## quantiles of the chi-squared distribution
    if(which == "all" || which == "qqchi2") {
        .qqplot(md, p, cutoff=cutoff, classic=TRUE, id.n=id.n, ...)
    }

    if(which == "all" || which == "tolEllipsePlot") {
        if(p == 2)
            .tolellipse(ccov = x, cutoff=cutoff, id.n=id.n, tol=tol, ...)
        else if(which != "all")
            warning("Warning: For tolerance ellipses the dimension must be 2!")
    }

    if(which == "tolEllipsePlot" || which == "pairs") {
        if(which == "tolEllipsePlot" & length(dim(data)) >= 2 && dim(data)[2] == 2){
            if(!is.null(rd)){
                .tolellipse(rcov=x, cutoff=cutoff, id.n=id.n, tol=tol, ...)
            }
        }else if(length(dim(data)) >= 2 && dim(data)[2] <= 10)
        {
            .rrpairs(x, ...)
        }else if(which != "all")
            warning("Warning: For tolerance ellipses the dimension must be less than 10!")
    }

    if(which == "all" || which == "screeplot") {
        myscreeplot(ccov=x)
    }
}) ## end { plot("CovClassic") }
