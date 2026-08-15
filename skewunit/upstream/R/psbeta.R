psbeta <-
function(q, delta=1, lower.tail=TRUE, log.p=FALSE)
{
if(delta<=0) stop("delta should be positive")
Fy<-pbeta(q, shape1=delta, shape2=delta, lower.tail=lower.tail, log.p=log.p)
}
