
dCBellBX<-function (x, a,lambda, log = FALSE)
{
	G=(1-exp(-x^2))^a
	g=2*a*x*exp(-x^2)*(1-exp(-x^2))^(a-1)
	pdf <- x
	pdf[log == FALSE] <- lambda*g*exp(lambda*G)*exp(exp(lambda*G)-1)/(exp(exp(lambda)-1)-1)
    pdf[log == TRUE] <- log(lambda)+log(g)+(lambda*G)+(exp(lambda*G)-1)-log((exp(exp(lambda)-1)-1))
	return(pdf)
}

