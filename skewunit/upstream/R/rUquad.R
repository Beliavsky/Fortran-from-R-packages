rUquad <-
function(n, a=0, b=1)
{
if (is.null(n)) stop("sample size must be specified")
if (round(n) != n | n <= 0) stop("sample size must be a positive integer")
alpha<-12/(b-a)^3
beta<-(a+b)/2
if(a>=b) stop("a must be less than b")
beta+cuberoot(3*runif(n)/alpha-(beta-a)^3)
}
