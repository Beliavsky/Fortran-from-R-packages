
#' @import stats

rCBellB<-function(n, a, b, k, lambda)
	{
    p <- runif(n, min = 0, max = 1)
    rn <- qCBellB(p, a, b, k, lambda)
    return(rn)
	}
