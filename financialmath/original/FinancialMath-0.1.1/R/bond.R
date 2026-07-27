bond=function(f,r,c,n,i,ic=1,cf=1,t=NA,plot=FALSE){

  all=list(f,r,c,n,i,ic,cf,t,plot)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(all,length) != 1)==T) stop("All inputs must be of length 1.")
  #NA
  nalist=list(f,r,c,n,i,ic,cf,plot)
  if(any(lapply(nalist,is.na)==T)) stop("Cannot input any variables except for t as NA.")
  #Numeric
  num2=list(f,r,c,n,i,ic,cf,t)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F)) stop("f, r, c, n, i, ic, cf, and t must be numeric.")
  #Logical
  stopifnot(is.logical(plot))

  NA.Neg=c(f,c,n,i,ic)
  #Positive
  if(any(NA.Neg<=0)) stop("All numeric variables must be positive.")
  if(!is.na(t) & (t<=0 | r<0)) stop("All numeric variables must be positive.")
  #Infinite
  if(any(NA.Neg==Inf)) stop("Cannot input any variables as infinite.")
  if(!is.na(t) & (t==Inf | r==Inf)) stop("Cannot input any variables as infinite.")

  if(n !=round(n)) stop("n must be a positive integer.")
  if(!is.na(t) & t>n) stop("t must be less than n.")

  eff.i=(1+i/ic)^(ic)-1
  int=(1+eff.i)^(1/cf)-1

  full=NA;clean=NA;p.d=NA;p=NA;prem=NA;dis=NA
  r=r/cf
  i=seq(0,2*int,by=.001)

  price=f*r*(1-(1+int)^(-n))/int+c*(1+int)^(-n)
  if(plot==T) price.g=f*r*(1-(1+i)^(-n))/i+c*(1+i)^(-n)
  if(!is.na(t) & t!=0) {
    p=f*r*(1-(1+int)^(-(n-t)))/int+c*(1+int)^(-(n-t))
    if(round(t) != t) {p=NA;d=floor(t);k=t-d
    full=(f*r*(1-(1+int)^(-(n-d)))/int+c*(1+int)^(-(n-d)))*(1+int)^k
    clean=full-k*f*r
    }
  }
  if(!is.na(t) & t!=0) p.d=abs(f*r-c*int)/(1+int)^(n-t+1)
  if(plot==T){
    plot(i,price.g,type="l",ylab="Price",xlab="Yield",main="Convexity of a Bond")
    abline(h=price,col="gray",lty=2)
    abline(v=int,col="gray",lty=2)  }
  n.v=seq(1,n,by=1)
  macd=(sum(f*r*n.v*(1+int)^(-n.v))+c*n*(1+int)^(-n))/price
  modd=macd/(1+int)
  macc=(sum(f*r*n.v*n.v/(1+int)^(n.v))+c*n*n*(1+int)^(-n))/price
  modc=(sum(f*r*n.v*(n.v+1)/(1+int)^(n.v+2))+c*n*(n+1)/(1+int)^(n+2))/price

  if(round(price,2)>round(c,2)) prem=price-c
  if(round(price,2)<round(c,2)) dis=c-price

  if(cf==1) {years=n;n=NA} else years=n/cf
  nom1=NA;nom2=NA
  coupon=f*r
  if(cf!=1) nom1=((1+eff.i)^(1/cf)-1)*cf
  if(ic!=1 & cf!=ic) nom2=((1+eff.i)^(1/ic)-1)*ic

  if(!is.na(t) & t==0) t=NA
  out=c(price,prem,dis,coupon,eff.i,nom1,nom2,n,years,macd,modd,macc,modc,t,p,full,clean,p.d)
  out=matrix(out,nrow=length(out))
  rownames(out)=c("Price","Premium","Discount","Coupon","Eff Rate",
                  paste("i^(",round(cf,2),")",sep=""),paste("i^(",round(ic,2),")",sep=""),
                  "Periods","Years","MAC D","MOD D","MAC C","MOD C",
                  paste("At Period ",round(t,2),":",sep=""),
                  "Price(t)","Full Price","Clean Price",
                  if(round(price,2)>round(c,2)) "Write-Down" else "Write-Up")
  na=apply(out, 1, function(x) all(is.na(x)))
  out=as.matrix(out[!na,])
  colnames(out)="Bond Summary"
  return(out)
}
