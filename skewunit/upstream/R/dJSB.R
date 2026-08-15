dJSB <-
function(x, delta=1, log=FALSE)
{
if(delta<=0) stop("delta should be positive")
eta<-log(x)-log1p(-x)
log.lik<-log(delta)-log(x)-log1p(-x)+dnorm(delta*eta, log=TRUE)
if(log == FALSE) fy <- exp(log.lik)
else fy <- log.lik
fy <- ifelse(x <= 0 | x>=1, 0, fy)
fy
}
