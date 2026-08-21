
#' @import stats

rCBellBX<-function(n, a,lambda)
	{
    p <- runif(n, min = 0, max = 1)
    rn <- qCBellBX(p, a, lambda)
    return(rn)
	}
