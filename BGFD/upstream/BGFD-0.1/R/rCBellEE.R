
#' @import stats

rCBellEE<-function(n, alpha,beta,lambda)
	{
    p <- runif(n, min = 0, max = 1)
    rn <- qCBellEE(p, alpha, beta, lambda)
    return(rn)
	}
