
#' @import stats

rBellEE<-function(n, alpha,beta,lambda)
	{
    p <- runif(n, min = 0, max = 1)
    rn <- qBellEE(p, alpha, beta, lambda)
    return(rn)
	}
