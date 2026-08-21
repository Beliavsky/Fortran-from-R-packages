
#' @import stats

rBellL<-function(n, b, q, lambda)
	{
    p <- runif(n, min = 0, max = 1)
    rn <- qBellL(p, b, q, lambda)
    return(rn)
	}
