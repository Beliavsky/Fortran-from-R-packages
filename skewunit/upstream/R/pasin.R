pasin <-
function(q, lower.tail=TRUE, log.p=FALSE)
{
Fy<-pbeta(q, shape1=1/2, shape2=1/2, lower.tail=lower.tail, log.p=log.p)
}
