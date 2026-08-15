dsbeta <-
function(x, delta=1, log=FALSE)
{
if(delta<=0) stop("delta should be positive")
log.lik<-dbeta(x, shape1=delta, shape2=delta, log=TRUE)
if(!log) fy <- exp(log.lik)
else fy <- log.lik
fy <- ifelse(x <= 0 | x>=1, 0, fy)
fy
}
