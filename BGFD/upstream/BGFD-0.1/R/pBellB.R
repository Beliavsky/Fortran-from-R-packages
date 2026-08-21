
pBellB<-function (x, a, b, k, lambda , log.p = FALSE, lower.tail = TRUE)
	{
	G=(1-((1+(x/a)^b))^(-k))
	cdf <- x
    cdf[log.p == FALSE & lower.tail == TRUE] <- (1-exp(-exp(lambda)*(1-exp(-lambda*G))))/(1-(exp(-exp(lambda)+1)))
	cdf[log.p == TRUE & lower.tail == TRUE] <- log(1-exp(-exp(lambda)*(1-exp(-lambda*G))))-log(1-(exp(-exp(lambda)+1)))
    cdf[log.p == FALSE & lower.tail == FALSE] <- ((exp(-exp(lambda)*(1-exp(-lambda*G))))-(exp(-exp(lambda)+1)))/(1-(exp(-exp(lambda)+1)))
    cdf[log.p == TRUE & lower.tail == FALSE] <- log((exp(-exp(lambda)*(1-exp(-lambda*G))))-(exp(-exp(lambda)+1)))-log(1-(exp(-exp(lambda)+1)))
    return(cdf)
	}
