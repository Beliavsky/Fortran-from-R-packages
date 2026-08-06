setMethod("getQuan", "PcaCov", function(obj) obj@n.obs)

##  The S3 version
PcaCov <- function (x, ...)
    UseMethod("PcaCov")

#' @export
PcaCov.formula <- function (formula, data = NULL, subset, na.action, ...)
{
    cl <- match.call()

    mt <- terms(formula, data = data)
    if (attr(mt, "response") > 0)
        stop("response not allowed in formula")
    mf <- match.call(expand.dots = FALSE)
    mf$... <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    ## this is not a 'standard' model-fitting function,
    ## so no need to consider contrasts or levels
    if (.check_vars_numeric(mf))
        stop("PCA applies only to numerical variables")

    na.act <- attr(mf, "na.action")
    mt <- attr(mf, "terms")
    attr(mt, "intercept") <- 0
    x <- model.matrix(mt, mf)

    res <- PcaCov.default(x, ...)

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("PcaCov")
    res@call <- cl

#    if (!is.null(na.act)) {
#        res$na.action <- na.act
#        if (!is.null(sc <- res$x))
#            res$x <- napredict(na.act, sc)
#    }

    res
}

#' @export
PcaCov.default <- function(x, k=ncol(x), kmax=ncol(x), cov.control = CovControlMcd(),
    scale=FALSE, signflip=TRUE, crit.pca.distances=0.975,
    trace=FALSE, ...)
{

    cl <- match.call()

    if(missing(x)){
        stop("You have to provide at least some data")
    }
    data <- as.matrix(x)
    n <- nrow(data)
    p <- ncol(data)

    if(n < p)
        stop("'PcaCov' can only be used with more units than variables")

    ## verify and set the input parameters: k and kmax
    myrank <- rankMM(x)
    kmax <- max(min(floor(kmax), myrank),1)
    if((k <- floor(k)) < 0)
        k <- 0
    else if(k > kmax) {
        warning(paste("The number of principal components k = ", k, " is larger then kmax = ", kmax, "; k is set to ", kmax,".", sep=""))
        k <- kmax
    }
######################################################################

    ##
    ## VT::05.08.2016
    ## If scale is TRUE/FALSE, this will be handled by .xpc()
    ##  otherwise call the undocumented function doScale() from robustbase -
    ##  it will behave properly and scale can be a vector or a function
    if(is.null(scale))
        scale <- vector('numeric', p) + 1
    else if(!is.logical(scale))
    {
        data <- doScale(data, center=NULL, scale=scale)
        scale <- data$scale
        data=data$x
    }

    ## VT::30.09.2009 - add the option for classic covariance estimates - if cov.control = NULL
    covx <- if(!is.null(cov.control)) restimate(cov.control, data) else Cov(data)
    covmat <- list(cov=getCov(covx), center=getCenter(covx), n.obs=covx@n.obs)

    ## VT::05.06.2016  - the call to princomp() replaced by an internal function
    ##  it will handle the case scale=TRUE and will return also the proper scores
    out <- .xpc(x, covmat=covmat, scale=scale, signflip=signflip)

## VT::11.28.2015: Choose the number of components k (if not specified)
##      (see mail of Klaus Nordhausen from 19.11.2015: the help says that the algorithm defines k)
##      before it was just k <- min(kmax, p), i.e. k=rank(X)
    if(k != 0)
        k <- min(k, p)
    else
    {
#        k <- min(kmax, p)
        ##
        ## Find the number of PC 'k'
        ## Use the test l_k/l_1 >= 10.E-3, i.e. the ratio of
        ## the k-th eigenvalue to the first eigenvalue (sorted decreasingly) is larger than
        ## 10.E/3 and the fraction of the cumulative dispersion is larger or equal 80%
        ##
        rk <- min(n, p)
        ev <- out$sdev^2
        test <- which(ev/ev[1] <= 1.E-3)
        k <- if(length(test) != 0)  min(min(rk, test[1]), kmax)
             else                   min(rk, kmax)

        cumulative <- cumsum(ev[1:k])/sum(ev)
        if(cumulative[k] > 0.8) {
            k <- which(cumulative >= 0.8)[1]
        }
        if(trace)
            cat("\n k, kmax, rank, p: ", k, kmax, rk, p, "\n")
        if(trace)
            cat("The number of principal components is defined by the algorithm. It is set to ", k,".\n", sep="")
    }

    center   <- getCenter(covx)
    scale    <- out$scale
    sdev     <- out$sdev
    loadings <- out$loadings[, 1:k, drop=FALSE]
    eigenvalues  <- (sdev^2)[1:k]
    scores   <- out$scores[, 1:k, drop=FALSE]
    eig0 <- sdev^2
    totvar0 <- sum(eig0)

######################################################################
    names(eigenvalues) <- NULL
    if(is.list(dimnames(x)))
    {
        rownames(scores) <- rownames(x)  # dimnames(scores)[[1]] <- dimnames(data)[[1]]
    }else
        dimnames(scores)[[1]] <- 1:n
    dimnames(scores)[[2]] <- paste("PC", seq_len(ncol(scores)), sep = "")
    dimnames(loadings) <- list(colnames(data), paste("PC", seq_len(ncol(loadings)), sep = ""))

    ## fix up call to refer to the generic, but leave arg name as 'formula'
    cl[[1]] <- as.name("PcaCov")
    res <- new("PcaCov", call=cl,
                            rank=myrank,
                            loadings=loadings,
                            eigenvalues=eigenvalues,
                            center=center,
                            scale=scale,
                            scores=scores,
                            k=k,
                            n.obs=n,
                            eig0=eig0,
                            totvar0=totvar0)

    ## Compute distances and flags
    res <- pca.distances(res, x, p, crit.pca.distances)
    return(res)
}

## A simplified version of princomp()
.xpc <- function (x, covmat, scale=FALSE, signflip=FALSE)
{
    ## x is always available and covmat is always a list
    ## scores is always TRUE (therefore we need x)
    ##
    if (any(is.na(match(c("cov", "n.obs"), names(covmat)))))
        stop("'covmat' is not a valid covariance list")
    cv <- covmat$cov
    n.obs <- covmat$n.obs
    cen <- covmat$center

    if(!is.numeric(cv))
        stop("PCA applies only to numerical variables")
    cor <- FALSE
    if(is.logical(scale))
    {
        if(scale)
        {
            cor <- TRUE
            sds <- sqrt(diag(cv))
            if(any(sds == 0))
                stop("cannot use 'cor = TRUE' with a constant variable")
            cv <- cv/(sds %o% sds)
        }
    }

    edc <- eigen(cv, symmetric = TRUE)
    ev <- edc$values
    if(any(neg <- ev < 0))
    {
        if(any(ev[neg] < -9 * .Machine$double.eps * ev[1L]))
            stop("covariance matrix is not non-negative definite")
        else
            ev[neg] <- 0
    }

    if(signflip)
        edc$vectors <- .signflip(edc$vectors)

    cn <- paste0("PC", 1L:ncol(cv))
    names(ev) <- cn
    dimnames(edc$vectors) <- list(dimnames(x)[[2L]], cn)
    sdev <- sqrt(ev)

    sc <- setNames(if(cor) sds
                   else if(!is.logical(scale)) scale
                   else rep.int(1, ncol(cv)), colnames(cv))
    scr <- as.matrix(doScale(x, center=cen, scale=sc)$x) %*% edc$vectors
    dimnames(scr) <- list(dimnames(x)[[1L]], cn)

    edc <- list(sdev=sdev, loadings=structure(edc$vectors, class="loadings"),
        center=cen, scale=sc, n.obs=n.obs, scores=scr)

    class(edc) <- "princomp"
    edc
}
