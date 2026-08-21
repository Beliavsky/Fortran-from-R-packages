
#' @import stats

rCBellE<-function(n, alpha,lambda)
	{
    p <- runif(n, min = 0, max = 1)
    rn <- qCBellE(p, alpha, lambda)
    return(rn)
	}
