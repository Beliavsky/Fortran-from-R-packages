
pCBellBX<-function (x,  a,lambda , log.p = FALSE, lower.tail = TRUE)
	{
	G=(1-exp(-x^2))^a
	cdf <- x
    cdf[log.p == FALSE & lower.tail == TRUE] <-  (exp(exp(lambda*G)-1)-1)/(exp(exp(lambda)-1)-1)        
	cdf[log.p == TRUE & lower.tail == TRUE]  <-  log(exp(exp(lambda*G)-1)-1)-log(exp(exp(lambda)-1)-1)  
	cdf[log.p == FALSE & lower.tail == FALSE] <- (exp(exp(lambda)-1)-exp(exp(lambda*G)-1))/((exp(exp(lambda)-1)-1))  
	cdf[log.p == TRUE & lower.tail == FALSE] <-  log(exp(exp(lambda)-1)-exp(exp(lambda*G)-1))-log((exp(exp(lambda)-1)-1))
    return(cdf)
	}
