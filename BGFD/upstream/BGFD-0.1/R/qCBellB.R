
qCBellB<-function(p, a, b, k, lambda,log.p = FALSE, lower.tail = TRUE){
	if (log.p == TRUE)
        p <- exp(p)
    if (lower.tail == FALSE)
        p <- 1 - p
	qf <- rep(NaN, length(p))
	t=(1/lambda*log(log(((p[p >= 0 & p <= 1])*(exp(exp(lambda)-1)-1))+1)+1))

	qf[p >= 0 & p <= 1] <- a*(((1-t)^(-1/k) -1)^(1/b))
    return(qf)
	}




