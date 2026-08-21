
qBellEW<-function(p, alpha, beta, theta, lambda, log.p = FALSE, lower.tail = TRUE){
	if (log.p == TRUE)
        p <- exp(p)
    if (lower.tail == FALSE)
        p <- 1 - p
	qf <- rep(NaN, length(p))
	t=(-1/lambda*log(1-((log(1-(p[p >= 0 & p <= 1])*(1-(exp(-exp(lambda)+1)))))/(-exp(lambda)))))

	qf[p >= 0 & p <= 1] <- (-1/alpha*log(1-(t)^(1/theta)))^(beta)
    return(qf)
	}




