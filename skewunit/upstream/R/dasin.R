dasin <-
function(x, log=FALSE)
{
log.lik<-dbeta(x, shape1=1/2, shape2=1/2, log=TRUE)
if(!log) fy <- exp(log.lik)
else fy <- log.lik
fy <- ifelse(x <= 0 | x>=1, 0, fy)
fy
}
