pskewunit<-function(x, lambda = 0, delta = 1, delta2 = 1, family1 = "asin", 
    family2 = "asin", log.p = FALSE)
{
	if (!any(family1 == c("asin", "sbeta", "Uquad", "triang", 
        "JSB"))) 
        stop("distribution is not recognized")
    if (!any(family2 == c("asin", "sbeta", "Uquad", "triang", 
        "JSB")) & !is.null(family2)) stop("distribution is not recognized")
    f <- get(paste("d", family1, sep = ""))
    F <- get(paste("p", family1, sep = ""))
    if(!is.null(family2)) G <- get(paste("p", family2, sep = ""))
    if (delta <= 0) 
        stop("delta should be positive")
    if (lambda < -1 | lambda > 1) 
        stop("lambda should be between -1 and 1")
	if(!is.null(family2))
	{
		p=c()
		for(i in 1:length(x))
		{
			p[i]=integrate(dskewunit, lambda = lambda, delta = delta, delta2 = delta2, family1 = family1, 
			    family2 = family2, lower=0, upper=x[i])$value
		}
	}
	if(is.null(family2))
	{
		p=F(x)
		if(any(family1 == c("sbeta", "JSB")))
		{
			p=F(x, delta=delta)
		}		
	}	
	p <- ifelse(x <= 0, 0, p)
	p <- ifelse(x >= 1, 1, p)
	if(log.p) p=log(p)
	p
}
