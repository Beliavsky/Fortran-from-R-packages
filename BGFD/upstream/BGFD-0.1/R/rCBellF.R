
#' @import stats

rCBellF<-function(n, a, b, lambda)
	{
    p <- runif(n, min = 0, max = 1)
    rn <- qCBellF(p, a, b, lambda)
    return(rn)
	}
