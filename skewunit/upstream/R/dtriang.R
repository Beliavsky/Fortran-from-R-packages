dtriang <-
function(x, log=FALSE)
{
log.lik<-ifelse(x<=1/2, log(4)+log(x), log(4)+log1p(-x))
if(log == FALSE) fy <- exp(log.lik)
else fy <- log.lik
fy <- ifelse(x <= 0 | x>=1, 0, fy)
fy
}
