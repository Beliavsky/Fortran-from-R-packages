pUquad <-
function(q, a=0, b=1, lower.tail=TRUE, log.p=FALSE)
{
alpha<-12/(b-a)^3
beta<-(a+b)/2
if(a>=b) stop("a must be less than b")
Fy<-log(alpha)-log(3)+log((q-beta)^3+(beta-a)^3)
Fy <- ifelse(q <= a, -Inf, Fy)
Fy <- ifelse(q >= b, 0, Fy)
if(log.p == FALSE) Fy <- exp(Fy)
if(lower.tail==FALSE) Fy<-1-Fy
Fy
}
