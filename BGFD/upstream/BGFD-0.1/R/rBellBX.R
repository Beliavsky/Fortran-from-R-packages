
#' @import stats

rBellBX<-function(n, a,lambda)
	{
    p <- runif(n, min = 0, max = 1)
    rn <- qBellBX(p, a, lambda)
    return(rn)
	}
