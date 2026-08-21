
#' @import stats

rBellEW<-function(n, alpha,beta,theta,lambda)
	{
    p <- runif(n, min = 0, max = 1)
    rn <- qBellEW(p, alpha, beta, theta, lambda)
    return(rn)
	}
