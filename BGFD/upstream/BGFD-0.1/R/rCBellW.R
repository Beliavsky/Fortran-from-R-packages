
#' @import stats

rCBellW<-function(n, alpha,beta,lambda)
	{
    p <- runif(n, min = 0, max = 1)
    rn <- qCBellW(p, alpha, beta, lambda)
    return(rn)
	}
