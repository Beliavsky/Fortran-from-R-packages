
#' @import stats

rCBellEW<-function(n, alpha,beta,theta,lambda)
	{
    p <- runif(n, min = 0, max = 1)
    rn <- qCBellEW(p, alpha, beta, theta, lambda)
    return(rn)
	}
