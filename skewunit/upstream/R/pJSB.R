pJSB <-
function(q, delta=1, lower.tail=TRUE, log.p=FALSE)
{
if(delta<=0) stop("delta should be positive")
eta<-log(q)-log1p(-q)
Fy<-pnorm(delta*eta, log.p=TRUE, lower.tail=TRUE)
Fy <- ifelse(q <= 0, -Inf, Fy)
Fy <- ifelse(q >= 1, 0, Fy)
if(log.p == FALSE) Fy <- exp(Fy)
if(lower.tail==FALSE) Fy<-1-Fy
Fy
}
