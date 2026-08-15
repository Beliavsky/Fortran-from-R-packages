ptriang <-
function(q, lower.tail=TRUE, log.p=FALSE)
{
Fy <- ifelse(q<=1/2, log(2)+2*log(q), log(2*q^2-(2*q-1)^2))
Fy <- ifelse(q <= 0, -Inf, Fy)
Fy <- ifelse(q >= 1, 0, Fy)
if(log.p == FALSE) Fy <- exp(Fy)
if(lower.tail==FALSE) Fy<-1-Fy
Fy
}
