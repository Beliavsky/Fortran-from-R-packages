
#' @import stats

rCBellL<-function(n, b, q, lambda)
	{
    p <- runif(n, min = 0, max = 1)
    rn <- qCBellL(p, b, q, lambda)
    return(rn)
	}
