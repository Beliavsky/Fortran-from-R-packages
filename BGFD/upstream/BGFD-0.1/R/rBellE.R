
#' @import stats

rBellE<-function(n, alpha,lambda)
	{
    p <- runif(n, min = 0, max = 1)
    rn <- qBellE(p, alpha, lambda)
    return(rn)
	}
