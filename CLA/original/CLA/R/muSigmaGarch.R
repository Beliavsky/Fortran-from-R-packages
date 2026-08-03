#### Yanhao Shi's implementation (in her MSc thesis) to get
#### (mu, Sigma)  from timeseries of assets.
#### ===========  Basically her  Index.ms2()  function
##
muSigmaGarch <- function(x,
                         formula = ~ garch(1,1),
                         cond.dist = "std", # <- standard t-distrib (w/ df = 5) ?
                         trace = FALSE, ...)
{
    if(length(d <- dim(x)) != 2) stop("'x' must be a numeric matrix or data frame (alike)")
    if(!requireNamespace("fGarch"))
        stop("muSigmaGarch() needs the 'fGarch' package installed")
    n <- d[1] # number of observations
    nA <- d[2]# number of assets
    if(any(d < 2)) stop("must have at least 2 observations of at least 2 assets")
    lr <- log(x[-1, ]/x[-n,]) # log return
    loss <- -lr
    garchFit <- fGarch :: garchFit
    gfit <- if(is.data.frame(loss))
                lapply(loss,    garchFit, formula=formula, cond.dist=cond.dist, trace=trace)
            else if(is.matrix(loss))
                apply(loss, 2L, garchFit, formula=formula, cond.dist=cond.dist, trace=trace)
            else stop("Invalid 'loss' (log return computation problem?)")
    gPRED <- fGarch::predict
    gpredict <- t(sapply(gfit, gPRED, n.ahead = 1))
    mu <- -unlist(gpredict[, "meanForecast"])
    sd <-  unlist(gpredict[, "standardDeviation"])
    x.std <- sapply(1:nA, function(i) lr[,i]/gfit[[i]]@sigma.t)
    ## diag.sd <- diag(sd)
    V <- cor(x.std, ...) # including names
    V[] <- sd * V * rep(sd, each = nA)
    ## = diag.sd %*% cor(x.std) %*% diag.sd
    list(mu = mu, covar = V)
}

