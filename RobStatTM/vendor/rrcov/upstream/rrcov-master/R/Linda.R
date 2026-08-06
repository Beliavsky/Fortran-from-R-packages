##  The S3 version
Linda <- function (x, ...) UseMethod("Linda")

#' @export
Linda.formula <- function(formula, data, ..., subset, na.action)
{
    m <- match.call(expand.dots = FALSE)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)
    xint <- match("(Intercept)", colnames(x), nomatch=0)
    if(xint > 0)
        x <- x[, -xint, drop=FALSE]
    res <- Linda.default(x, grouping, ...)

##    res$terms <- Terms

    ## fix up call to refer to the generic, but leave arg name as 'formula'
    cl <- match.call()
    cl[[1]] <- as.name("Linda")
    res@call <- cl

##    res$contrasts <- attr(x, "contrasts")
##    res$xlevels <- .getXlevels(Terms, m)
##    res$na.action <- attr(m, "na.action")

    res
}

#' @export
Linda.default <- function(x,
                 grouping,
                 prior = proportions,
                 tol = 1.0e-4,
                 method = c("mcd", "mcdA", "mcdB", "mcdC", "fsa", "mrcd", "ogk"),
                 alpha=0.5,
                 l1med=FALSE,
                 cov.control,
                 trace=FALSE, ...)
{
    if(is.null(dim(x)))
        stop("x is not a matrix")

    method <- match.arg(method)
    xcall <- match.call()
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)

    if(length(grouping) == 1) {
        ## this is the number of groups and the groups are of equal size
        ng = grouping
        ni = n/ng
        if(ng*ni < n)
            stop("nrow(x) is not divisible by the number of groups")
        grouping <- rep(0,0)
        for(i in 1:ng)
            grouping <- c(grouping, rep(i,ni))
    }else if(length(grouping) > 1 && length(grouping) < n) {
        ## grouping contains a vector with the group sizes
        ng <- length(grouping)
        if(sum(grouping) != n)
            stop("nrow(x) is not equal to n1+n2+...+nn")

        gx <- rep(0,0)
        for(i in 1:ng)
            gx <- c(gx, rep(i,grouping[i]))
        grouping <- gx
    }

    if(n != length(grouping))
        stop("nrow(x) and length(grouping) are different")

    g <- as.factor(grouping)
    lev <- lev1 <- levels(g)
    counts <- as.vector(table(g))

    if(!missing(prior)) {
        if(any(prior < 0) || round(sum(prior), 5) != 1)
            stop("invalid prior")
        if(length(prior) != nlevels(g))
            stop("prior is of incorrect length")
        prior <- prior[counts > 0]
    }
    if(any(counts == 0)) {
        warning(paste("group(s)", paste(lev[counts == 0], collapse=" "),"are empty"))
        lev1 <- lev[counts > 0]
        g <- factor(g, levels=lev1)
        counts <- as.vector(table(g))
    }
    proportions <- counts/n
    ng <- length(proportions)
    names(g) <- NULL
    names(prior) <- levels(g)

    if(missing(cov.control))
        cov.control <- if(method == "mrcd") CovControlMrcd(alpha=alpha)
                       else if(method == "ogk") CovControlOgk()
                       else if(method %in% c("mcd", "mcdA", "mcdB", "mcdC")) CovControlMcd(alpha=alpha)
                       else NULL

    if(method == "mcdC" && !inherits(cov.control, "CovControlMcd"))
        stop("Method 'C' is defined only for MCD estimates!")

    if(method == "mcdA" && !inherits(cov.control, "CovControlMcd"))
        stop("Method 'A' is defined only for MCD estimates!")

    if(method == "mrcd" && !inherits(cov.control, "CovControlMrcd"))
        stop("Method 'mrcd' is defined only for MRCD estimates!")

    if(method == "ogk" && !inherits(cov.control, "CovControlOgk"))
        stop("Method 'ogk' is defined only for OGK estimates!")

    if(inherits(cov.control, "CovControlMrcd"))
        method <- "mrcd"
    if(inherits(cov.control, "CovControlOgk"))
        method <- "ogk"

    if(method == "fsa"){
        if(nrow(x) > 5000 | ncol(x) > 100)
            stop("Method 'fsa' can handle at most 5000 cases and 100 variables!")
        xcov <- .wcovMwcd(x, grouping, alpha=alpha, trace=trace)
    } else if(method == "mcdA"){
        xcov <- .wcovMcd(x, grouping, method="A", cov.control=cov.control)
    } else if(method == "mcd" || method == "mrcd" || method == "mcdB" || method == "ogk"){
        xcov <- .wcovMcd(x, grouping, method="B", l1med=l1med, cov.control=cov.control)
    } else if(method == "mcdC"){
        xcov <- .wcovMcd(x, grouping, method="C", l1med=l1med, cov.control=cov.control)
    } else {
        stop(paste("Unknown method called: ", method))
    }

    ##  VT::27.11.2019
    ##  inv <- solve(xcov$wcov)
    inv <- if(!is.null(xcov$winv)) xcov$winv
           else if(!.isSingular(xcov$wcov))  solve(xcov$wcov)
           else .pinv(xcov$wcov)

    ldf <- xcov$means %*% inv
    ldfconst <- diag(log(prior) - ldf %*% t(xcov$means)/2)
    if(!is.null(xcov$raw.means) && !is.null(xcov$raw.wcov)){
        raw.means <- xcov$raw.means
        raw.wcov <- xcov$raw.wcov
        ## inv <- solve(raw.wcov)
        inv <- if(!.isSingular(raw.wcov))  solve(raw.wcov) else .pinv(raw.wcov)
        raw.ldf <- raw.means %*% inv
        raw.ldfconst <- diag(log(prior) - raw.ldf %*% t(raw.means)/2)
    }else{
        raw.means <- xcov$means
        raw.wcov <- xcov$wcov
        raw.ldf <- ldf
        raw.ldfconst <- ldfconst
    }
    return (new("Linda", call=xcall, prior=prior, counts=counts,
                 center=xcov$means, cov=xcov$wcov, ldf = ldf, ldfconst = ldfconst,
                 method=method, l1med=l1med, X=x, grp=g,
                 covobj=xcov$covobj, control=cov.control))
}


.wcovMcd <- function(x, grouping, method = c("A", "B", "C"), alpha=0.5, l1med=FALSE, cov.control){
    xcall <- match.call()
    method <- match.arg(method)
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    dimn <- dimnames(x)

    if(!is.factor(g <- grouping))
        g <- as.factor(grouping)
    lev <- levels(g)
    counts <- as.vector(table(g))
    if(any(counts == 0)) {
        stop(paste("group(s)", paste(lev[counts == 0], collapse=" "),"are empty"))
    }
    ng <- length(counts/n)

    # compute group means and covariance matrices for each group
    mX <- matrix(0,ng,p)
    covX <- array(0,c(p,p,ng))
    sumwt <- array(0,ng)
    raw.mX <- matrix(0,ng,p)
    raw.covX <- array(0,c(p,p,ng))
    raw.sumwt <- array(0,ng)

    if(method == "A" | !l1med)
    {
        for(i in 1:ng){
            mcd <- if(!is.null(cov.control)) restimate(cov.control, x[which(g == lev[i]),]) else Cov(x[which(g == lev[i]),])
            mX[i,] <- getCenter(mcd)
            if(inherits(mcd, "CovMcd"))
            {
                sumwt[i] <- sum(mcd@wt)
                covX[,,i] <- getCov(mcd) * sumwt[i]
                raw.mX[i,] <- mcd@raw.center
                raw.sumwt[i] <- length(mcd@best)
                raw.covX[,,i] <- mcd@raw.cov * raw.sumwt[i]
            }
        }
    } else
    {
        for(i in 1:ng){
            mX[i,] <- l1median(x[which(g == lev[i]),])
        }
    }
    if(method == "A"){
        #Method A: pool the covariance matrices
        wcov <- matrix(0,p,p)
        for(i in 1:ng){
            wcov <- wcov + covX[,,i]
        }
        wcov <- wcov/(sum(sumwt)-ng)
        final.xcov <- list(wcov=wcov, means=mX, covobj=NULL)

        ## pool the raw estimates
        wcov <- matrix(0,p,p)
        for(i in 1:ng){
            wcov <- wcov + raw.covX[,,i]
        }
        wcov <- wcov/(sum(raw.sumwt)-ng)
        mX <- raw.mX
        mah <- mahalanobis(x-mX[g,], rep(0,p), wcov)
        weights <- ifelse(mah< qchisq(0.975, p), 1, 0)
        method <- "mcd-A"
    }else if(method == "B"){
        #Method B: center the data and compute the covariance matrix
        #of the centered data

        mcd <- if(!is.null(cov.control)) restimate(cov.control, x - mX[g,]) else Cov(x - mX[g,])
        winv <- if("icov" %in% slotNames(mcd)) mcd@icov else NULL
        final.xcov <- list(wcov=getCov(mcd), winv=winv, means=t(t(mX)+getCenter(mcd)), covobj=mcd)

        if(inherits(mcd, "CovMcd"))
        {
            wcov <- mcd@raw.cov
            mX <- t(t(mX)+mcd@raw.center)
            mah <- mcd@raw.mah
            weights <- mcd@raw.wt
        }else
        {
            wcov <- getCov(mcd)
            mX <- t(t(mX)+getCenter(mcd))
            mah <- weights <- NULL
        }

        method <- "mcd-B"
    }else if(method == "C"){

        ## Method C is more or less the same as Method B, but the means
        ##  are computed as the means of the Hi observations from mcd$best
        ##  and the covariance matrix is the raw.cov
        ##
        ## Center the data and compute the means and covariance matrix
        ## of the centered data as in Method B
        mcd <- restimate(cov.control, x - mX[g,])

        ## compute the group means as the means of these observations which are in
        ##  mcd@best, i.e. in case of two groups, if H=best
        ##  partition H in H1 and H2 and compute mi as the mean of the Hi observations
        ##  respectively.
        ##  Take the raw covariance matrix as within group cov
        for(i in 1:ng){
            tmpc <- cov.wt(as.matrix(x[mcd@best[which(mcd@best%in% which(g == lev[i]))],]))
            mX[i,] <- tmpc$center
        }
        wcov <- mcd@raw.cov
        mah <- mahalanobis(x-mX[g,], rep(0,p), wcov)
        weights <- ifelse(mah< qchisq(0.975, p), 1, 0)
        final.xcov <- .wcov.wt(x, g, weights)

        method <- "mcd-C"
    }else{
        stop("Unknown method specified: ", method)
    }

    dimnames(wcov) <- list(dimn[[2]], dimn[[2]])
    dimnames(mX) <- list(levels(g), dimn[[2]])
    dimnames(final.xcov$wcov) <- list(dimn[[2]], dimn[[2]])
    dimnames(final.xcov$means) <- list(levels(g), dimn[[2]])

    ans <- list(call=xcall, means=final.xcov$means, wcov=final.xcov$wcov, winv=final.xcov$winv,
                method=method,
                raw.means=mX, raw.wcov=wcov, raw.mah=mah, raw.wt=weights,
                covobj=final.xcov$covobj)

    class(ans) <- "wcov"
    return(ans)
}

.wcovMwcd <- function(x, grouping, alpha=0.5, trace=0){

    quan.f <- function(alpha, n, rk)
    {
        quan <- floor(2 * floor((n+rk+1)/2) - n + 2 * (n - floor((n+rk+1)/2)) * alpha)
        return(quan)
    }

    xcall <- match.call()
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    quan <- quan.f(alpha, n, p)
    dimn <- dimnames(x)

    if(!is.factor(g <- grouping))
        g <- as.factor(grouping)
    lev <- levels(g)
    counts <- as.vector(table(g))
    if(any(counts == 0)) {
        stop(paste("group(s)", paste(lev[counts == 0], collapse=" "),"are empty"))
    }
    ng <- length(counts/n)
    grouping <- as.integer(g)

    # transpose x, as necessary for the FSADA subroutine
    X <- t(x)

    storage.mode(X) <- "double"
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"

    storage.mode(ng) <- "integer"                   # NG
    storage.mode(grouping) <- "integer"             # IGRP

    xmean <- matrix(0, nrow = p * ng, ncol = 1)     # XM
    storage.mode(xmean) <- "double"
    xcov <- matrix(0, nrow = p * p, ncol = 1)       # XC
    storage.mode(xcov) <- "double"

    storage.mode(counts) <- "double"                # XND

    storage.mode(quan) <- "integer"                 # IHALF
    nsamp <- 0
    storage.mode(nsamp) <- "integer"                # NIT

    inbest <- matrix(10000, nrow = quan, ncol = 1)  # IDSAV
    storage.mode(inbest) <- "integer"
    seed <- 0
    storage.mode(seed) <- "integer"                 # MRAND
    deter <- 0                                      # DETC
    storage.mode(deter) <- "double"
    ierr <- 0                                       # IERR
    storage.mode(ierr) <- "integer"


    mwcd <- .Fortran("fsada",
        X,
        n,
        p,
        ng,
        grouping,
        xmean = xmean,
        xcov = xcov,
        counts,
        quan,
        nsamp,
        best = inbest,
        seed,
        ierr = ierr,
        detc = deter,
        itrace=as.integer(trace),
        PACKAGE="rrcov")

    xcov <- mwcd$xcov
    xmean <- mwcd$xmean
    dim(xcov) <- c(p, p)
    dim(xmean) <- c(p,ng)
    xmean <- t(xmean)
    dimnames(xcov) <- list(dimn[[2]], dimn[[2]])
    dimnames(xmean) <- list(levels(g), dimn[[2]])

    ## Compute the consistency correction factor for the raw MCD
    ## (see calfa in croux and haesbroeck)
    qalpha <- qchisq(quan/n, p)
    calphainvers <- pgamma(qalpha/2, p/2 + 1)/(quan/n)
    calpha <- 1/calphainvers
    correct <- 1

    if(p == 1)
    {
        scale <- sqrt(calpha) * as.double(xcov) * sqrt(correct)
    }else
    {
        ## Apply correction factor to the raw estimates and use them to compute weights
        xcov <- calpha * xcov * correct
        mah <- mahalanobis(x-xmean[g,], rep(0,p), xcov)
        weights <- ifelse(mah< qchisq(0.975, p), 1, 0)
        final.cov <- .wcov.wt(x, g, weights)
    }

    ans <- list(call=xcall,
                raw.means=xmean, raw.wcov=xcov, raw.mah=mah, raw.wt=weights,
                means=final.cov$means, wcov=final.cov$wcov, covobj=NULL)
    class(ans) <- "wcov"
    return(ans)
}
