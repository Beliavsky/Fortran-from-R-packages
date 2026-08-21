
#' @import stats

rBellF<-function(n, a, b, lambda)
	{
    p <- runif(n, min = 0, max = 1)
    rn <- qBellF(p, a, b, lambda)
    return(rn)
	}
