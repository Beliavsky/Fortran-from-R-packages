setClass("PcaLocantore", representation(delta = "numeric",
                                    quan = "numeric"),
                                contains="PcaRobust")

setMethod("getQuan", "PcaLocantore", function(obj) obj@n.obs)

##  The S3 version
PcaLocantore <- function (x, ...)
    UseMethod("PcaLocantore")

#' @export
PcaLocantore.formula <- function (formula, data = NULL, subset, na.action, ...)
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

    res <- PcaLocantore.default(x, ...)

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("PcaLocantore")
    res@call <- cl

#    if (!is.null(na.act)) {
#        res$na.action <- na.act
#        if (!is.null(sc <- res$x))
#            res$x <- napredict(na.act, sc)
#    }

    res
}

#' @export
PcaLocantore.default <- function(x, k=ncol(x), kmax=ncol(x), delta = 0.001, na.action = na.fail,
    scale=FALSE, signflip=TRUE, crit.pca.distances=0.975, trace=FALSE, ...)
{

    cl <- match.call()

    if(missing(x)){
        stop("You have to provide at least some data")
    }
    y <- data <- x <- as.matrix(x)
    n <- nrow(y)
    p <- ncol(y)

    ##
    ## verify and set the input parameters: k and kmax
    ##
    myrank <- rankMM(x)
    kmax <- max(min(floor(kmax), myrank),1)
    if((k <- floor(k)) < 0)
        k <- 0
    else if(k > kmax) {
        warning(paste("The number of principal components k = ", k, " is larger then kmax = ", kmax, "; k is set to ", kmax,".", sep=""))
        k <- kmax
    }

    if(k != 0)
        k <- min(k, ncol(data))

######################################################################

    ## VT::4.08.2016
    ##  if scale = TRUE, scale with mad, otherwise scale can be a vector
    ##      of length p, or a function
    ##  use the undocumented function doScale() from robustbase
    ##
    if(is.logical(scale))
        scale <- if(scale) mad else NULL
    if(!is.null(scale))
    {
        x.scaled <- doScale(x, center=NULL, scale=scale)
        sc <- x.scaled$scale
        data <- x.scaled$x
    }else
    {
        sc <- vector('numeric', p) + 1
        data <- x
    }

##    if(is.logical(scale))
##    {
##        if(scale)
##        {
##            sc <- apply(x, 2L, sd)
##            data <- sweep(data, 2L, sc, "/", check.margin=FALSE)
##        }
##    }
##    else if(is.numeric(scale) && length(scale) == p)
##    {
##        data <- sweep(data, 2L, scale, "/", check.margin = FALSE)
##        sc <- scale
##    }
##    else stop("length of 'scale' must equal the number of columns of 'x'")
##    attr(data, "scaled:scale") <- sc

    ## Center with the spatial median and project on the sphere
    spa = spatial.median(data, delta)
    mu = spa$mu
    ep = spa$ep
    tt = matrix(mu, n, p, byrow=TRUE)
    data = data-tt
    for(i in 1:n)
    {
        z = sqrt(sum((data[i,  ])^2))
        y[i,  ] = 0 * data[i,  ]
        if(z > ep)
        {
            y[i,  ] = (data[i,  ]  )/z
        }
    }

    ## no scaling - we have already scaled with MAD
    if(k == 0)      # let PcaClassic guess the number of components to select
    {
        out = PcaClassic(y, k=k, kmax=kmax, scale=FALSE, signflip=signflip, ...)
        k <- out@k
    }

    ## This will return all components, later we will select the first k
    out = PcaClassic(y, scale=FALSE, signflip=signflip, ...)

    scores = y %*% out@loadings             # these are (slightly) diferent from the scores returned by PcaClassic
                                            #   because PcaClassic will center by the mean the already centered data
    # use these scores to compute the standard deviations (mad)
    sdev = apply(scores, 2, "mad")
    names2 = names(sdev)
    orsdev = order(sdev)            # sort the sdevs (although almost always they will be sorted)
    orsdev = rev(orsdev)            # use them to sort the loadings, etc.
    sdev  = sdev[orsdev]
    scores  = scores[,orsdev, drop=FALSE]
    loadings = out@loadings[,orsdev, drop=FALSE]

    names(sdev)=names2
    dimnames(scores)[[2]]=names2
    dimnames(loadings)[[2]]=names2

    ## scale is a vector of 1s with length p (if FALSE), a vector with column mad (if TRUE) or
    ##      any other vector with length p (posibly computed by doScale())
    scale       <- sc
    center      <- doScale(t(as.vector(mu)), center=NULL, scale=1/scale)$x      # rescale back to the original data
    scores      <- doScale(x, center, scale)$x %*% loadings                     # compute the scores
    scores      <- scores[, 1:k, drop=FALSE]                                    # select only fist k
    loadings    <- loadings[, 1:k, drop=FALSE]
    eigenvalues <- (sdev^2)[1:k]
    eig0 <- sdev^2
    totvar0 <- sum(eig0)

######################################################################
    names(eigenvalues) <- NULL
    if(is.list(dimnames(data)))
    {
        ##dimnames(scores)[[1]] <- dimnames(data)[[1]]
        rownames(scores) <- rownames(data)
    }
    dimnames(scores)[[2]] <- as.list(paste("PC", seq_len(ncol(scores)), sep = ""))
    dimnames(loadings) <- list(colnames(data), paste("PC", seq_len(ncol(loadings)), sep = ""))

    ## fix up call to refer to the generic, but leave arg name as <formula>
    cl[[1]] <- as.name("PcaLocantore")
    res <- new("PcaLocantore", call=cl,
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

## computes the spatial median
spatial.median <- function(x, delta)
{
    dime = dim(x)
    n=dime[1]
    p=dime[2]
    delta1=delta*sqrt(p)
    mu0=apply(x,2,median)
    h=delta1+1
    tt=0
    while(h>delta1)
    {
        tt=tt+1
        TT=matrix(mu0,n,p,byrow=TRUE)
        U=(x-TT)^2
        w=sqrt(apply(U,1,sum))
        w0=median(w)
        ep=delta*w0

        z=(w<=ep)
        w[z]=ep
        w[!z]=1/w[!z]
        w=w/sum(w)
        x1=x
        for(i in 1:n)
            x1[i,]=w[i]*x[i,]
        mu=apply(x1,2,sum)
        h=sqrt(sum((mu-mu0)^2))
        mu0=mu
    }
    out=list(mu=mu0,ep=ep)
    out
}
