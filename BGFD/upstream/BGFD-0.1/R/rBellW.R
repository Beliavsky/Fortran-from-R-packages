
#' @import stats

rBellW<-function(n, alpha,beta,lambda)
	{
    p <- runif(n, min = 0, max = 1)
    rn <- qBellW(p, alpha, beta, lambda)
    return(rn)
	}
