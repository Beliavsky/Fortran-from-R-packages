dUquad <-
function(x, a=0, b=1, log=FALSE)
{
alpha<-12/(b-a)^3
beta<-(a+b)/2
if(a>=b) stop("a must be less than b")
log.lik<-log(alpha)+2*log(abs(x-beta))
if(log == FALSE) fy <- exp(log.lik)
else fy <- log.lik
fy <- ifelse(x <= a | x>=b, 0, fy)
fy
}
