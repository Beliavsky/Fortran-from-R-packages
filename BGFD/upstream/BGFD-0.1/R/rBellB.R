
#' @import stats

rBellB<-function(n, a, b, k, lambda)
	{
    p <- runif(n, min = 0, max = 1)
    rn <- qBellB(p, a, b, k, lambda)
    return(rn)
	}
