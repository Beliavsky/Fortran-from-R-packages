rskewunit <-
function(n, lambda=0, delta=1, delta2=1, family1="asin", family2="asin")
{
 if (!any(family1 == c("asin", "sbeta", "Uquad", "triang", 
        "JSB"))) 
        stop("distribution is not recognized")
    if (!any(family2 == c("asin", "sbeta", "Uquad", "triang", 
        "JSB")) & !is.null(family2)) 
        stop("distribution is not recognized")
r<-get(paste("r",family1, sep=""))
x <- r(n)
if(family1 == "JSB" | family1 =="sbeta") x<-r(n, delta)
if(!is.null(family2))
{
G<-get(paste("p",family2, sep=""))
if (is.null(n)) stop("sample size must be specified")
if (round(n) != n | n <= 0) stop("sample size must be a positive integer")
if(delta<=0) stop("delta should be positive")
if(delta2<=0) stop("delta2 should be positive")
x=c();i=1
while(i<=n)
{
y <- r(1)
if(family1 == "JSB" | family1 =="sbeta") y<-r(1, delta)
prob<-G(lambda*(y-0.5)+0.5)
if(family2 == "JSB" | family1 =="sbeta") prob<-G(lambda*(y-0.5)+0.5, delta2)
if(runif(1)<=prob/2){x[i]=y; i=i+1}
}
}
x
}
