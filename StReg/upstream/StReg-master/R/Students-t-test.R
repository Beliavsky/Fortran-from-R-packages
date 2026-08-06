Student <- function(res,Data,Ct,s2,v,lag,B1)
{

if (ncol(res)>1) {p <- ncol(Data)}
if (ncol(res)==1) {p <- 1}
v <- ifelse(v<=1,1,v)
ad <- matrix(nrow=p,ncol=2)
for(k in 1:p)
{
u.hat <- res[,k]/sqrt(((v/(v+ncol(B1)-2))*(s2[k,k]))*Ct) 
AD <- ad.test(u.hat,pt,df=v+ncol(B1))
ad[k,] <- cbind(AD[1]$statistic, AD[2]$p.value)
}
colnames(ad) <- c("AD","p-value")

return(list(ad=ad))

}
