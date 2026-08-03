##setGeneric("PcaClassic", function(x, ...) standardGeneric("PcaClassic"))
##setMethod("PcaClassic", "formula", PcaClassic.formula)
##setMethod("PcaClassic", "ANY", PcaClassic.default)

setMethod("getQuan", "PcaClassic", function(obj) obj@n.obs)

##  The S3 version
PcaClassic <- function (x, ...) UseMethod("PcaClassic")

#' @export
PcaClassic.formula <- function (formula, data = NULL, subset, na.action, ...)
{
    cl <- match.call()

    mt <- terms(formula, data = data)
    if (attr(mt, "response") > 0)
        stop("response not allowed in formula")
    mf <- match.call(expand.dots = FALSE)
    mf$... <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    ## this is not a `standard' model-fitting function,
    ## so no need to consider contrasts or levels
    if (.check_vars_numeric(mf))
        stop("PCA applies only to numerical variables")

    na.act <- attr(mf, "na.action")
    mt <- attr(mf, "terms")
    attr(mt, "intercept") <- 0
    x <- model.matrix(mt, mf)

    res <- PcaClassic.default(x, ...)

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("PcaClassic")
    res@call <- cl

#    if (!is.null(na.act)) {
#        res$na.action <- na.act
#        if (!is.null(sc <- res$x))
#            res$x <- napredict(na.act, sc)
#    }

    res
}

#' @export
PcaClassic.default <- function(x, k=ncol(x), kmax=ncol(x),
    scale=FALSE, signflip=TRUE, crit.pca.distances=0.975, trace=FALSE, ...)
{
    cl <- match.call()

    if(missing(x)){
        stop("You have to provide at least some data")
    }
    data <- as.matrix(x)
    n <- nrow(data)
    p <- ncol(data)

    Xsvd <- .classPC(data, scale=scale, signflip=signflip, scores=TRUE)
    if(Xsvd$rank == 0) {
        stop("All data points collapse!")
    }

    myrank <- Xsvd$rank

    if(is.logical(scale) && !scale)         # no scaling, the defult
        Xsvd$scale <- vector('numeric', p) + 1

    if(trace)
    {
        cat("\nDimension of the input matrix x:\n", dim(x))
        cat("\nInput parameters [k, kmax, rank(x)]: ", k, kmax, Xsvd$rank, "\n")
    }

    ##
    ## verify and set the input parameters: k and kmax
    ##
    kmax <- max(min(kmax, Xsvd$rank),1)

    if((k <- floor(k)) < 0)
        k <- 0
    else if(k > kmax) {
        warning(paste("The number of principal components k = ", k, " is larger then kmax = ", kmax, "; k is set to ", kmax,".", sep=""))
        k <- kmax
    }

    if(k != 0)
        k <- min(k, ncol(data))
    else {
##        k <- min(kmax, ncol(data))
##        if(trace)
##            cat("The number of principal components is set to ", k, ".\n", sep="")

        ## Find the number of PC 'k'
        ## Use the test l_k/l_1 >= 10.E-3, i.e. the ratio of
        ## the k-th eigenvalue to the first eigenvalue (sorted decreasingly) is larger than
        ## 10.E/3 and the fraction of the cumulative dispersion is larger or equal 80%
        ##
        test <- which(Xsvd$eigenvalues/Xsvd$eigenvalues[1] <= 1.E-3)
        k <- if(length(test) != 0)  min(min(Xsvd$rank, test[1]), kmax)
             else                   min(Xsvd$rank, kmax)

        cumulative <- cumsum(Xsvd$eigenvalues[1:k])/sum(Xsvd$eigenvalues)
        if(cumulative[k] > 0.8) {
            k <- which(cumulative >= 0.8)[1]
        }

        if(trace)
            cat("\n k, kmax, rank, p: ", k, kmax, Xsvd$rank, ncol(data), "\n")
        if(trace)
            cat("The number of principal components is defined by the algorithm. It is set to ", k,".\n", sep="")

    }

    if(trace)
        cat("\nTo be used [k, kmax, ncol(data), rank(data)]=",k, kmax, ncol(data), Xsvd$rank, "\n")

    loadings    <- Xsvd$loadings[, 1:k, drop=FALSE]
    eigenvalues <- as.vector(Xsvd$eigenvalues[1:k])
    center      <- as.vector(Xsvd$center)
    scores      <- Xsvd$scores[, 1:k, drop=FALSE]
    scale       <- Xsvd$scale
    eig0 <- as.vector(Xsvd$eigenvalues)
    totvar0 <- sum(eig0)

    if(is.list(dimnames(data)) && !is.null(dimnames(data)[[1]]))
    {
        dimnames(scores)[[1]] <- dimnames(data)[[1]]
    } else {
        dimnames(scores)[[1]] <- 1:n
    }
    dimnames(scores)[[2]] <- as.list(paste("PC", seq_len(ncol(scores)), sep = ""))
    dimnames(loadings) <- list(colnames(data), paste("PC", seq_len(ncol(loadings)), sep = ""))

    ## fix up call to refer to the generic, but leave arg name as 'formula'
    cl[[1]] <- as.name("PcaClassic")
    res <- new("PcaClassic", call=cl,
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
    res <- pca.distances(res, data, Xsvd$rank, crit.pca.distances)
    return(res)
}
