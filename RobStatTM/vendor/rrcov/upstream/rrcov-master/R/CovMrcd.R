CovMrcd <- function(x,
                   alpha=control@alpha,
                   h=control@h,
                   maxcsteps=control@maxcsteps,
                   initHsets=NULL, save.hsets=FALSE,
                   rho=control@rho,
                   target=control@target,
                   maxcond=control@maxcond,
                   trace=control@trace,
                   control=CovControlMrcd())
{
    if(is.data.frame(x))
        x <- data.matrix(x, rownames.force=FALSE)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
            dimnames = list(names(x), deparse(substitute(x))))

    ## drop all rows with missing values (!!) :
    ok <- is.finite(x %*% rep.int(1, ncol(x)))
    x <- x[ok, , drop = FALSE]
    if(!length(dx <- dim(x)))
        stop("All observations have missing values!")
    n <- dx[1]; p <- dx[2]
    dimn <- dimnames(x)

    ## VT::18.07.2022 - maxcond passed to the lowlevel function
    mcd <- .detmrcd (x, alpha=alpha, h=h, hsets.init = initHsets,
		      save.hsets=save.hsets, # full.h=full.h,
		      rho=rho, maxcond=maxcond, target=if(target=="identity") 0 else 1,
              maxcsteps=maxcsteps,
              trace=as.integer(trace))

    alpha <- mcd$alpha
    h <- mcd$h
    ans <- list(call = match.call(), method = sprintf("MRCD(alpha=%g ==> h=%d)", alpha, h))
	ans$method <- paste("Minimum Regularized Covariance Determinant", ans$method)

    ans$cov <- mcd$initcovariance
    ans$center <- as.vector(mcd$initmean)

    ans$n.obs <- n
    ans$best <- sort(as.vector(mcd$best))
    ans$alpha <- alpha
    ans$quan <- h
    ans$crit <- mcd$mcdestimate
    ans$mah <- mahalanobis(x, mcd$initmean, mcd$icov, inverted=TRUE)

	if(length(dimn[[1]]))
	    dimnames(x)[[1]] <- dimn[[1]][ok]
	else
	    dimnames(x) <- list(seq(along = ok)[ok], NULL)

    ans$X <- x
    if(trace)
        cat(ans$method, "\n")

    ans <- c(ans, mcd[c("calpha", "iBest","n.csteps", if(save.hsets) "initHsets", "icov","rho", "target")])
    class(ans) <- "mcd"

    if(!is.null(nms <- dimn[[2]])) {
        dimnames(ans$cov) <- list(nms, nms)
        dimnames(ans$icov) <- list(nms, nms)
        names(ans$center) <- nms
    }

    new("CovMrcd",
        call= ans$call,
        crit=ans$crit,
        cov=ans$cov,
        icov=ans$icov,
        rho=ans$rho,
        target=ans$target,
        center=ans$center,
        n.obs=ans$n.obs,
        mah = ans$mah,
        X = ans$X,
        method=ans$method,
        best=ans$best,
        alpha=ans$alpha,
        quan=ans$quan,
        cnp2 = ans$calpha)
}
