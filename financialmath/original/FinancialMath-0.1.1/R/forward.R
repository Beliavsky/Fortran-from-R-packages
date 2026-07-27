forward=function(S,t,r,position,div.structure="none",dividend=NA,df=1,D=NA,k=NA,plot=FALSE){

  all=list(S,t,r,position,div.structure,dividend,df,D,k,plot)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(all,length) != 1)==T) stop("All inputs must be of length 1.")
  #Numeric
  num2=list(S,t,r,dividend,df,D,k)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F))stop("S, t, r, dividend, df, D, and k must be numeric.")
  #NA
  nalist=list(S,t,r,position,div.structure,df,plot)
  if(any(lapply(nalist,is.na)==T)) stop("S, t, r, position, div.structure, df, and plot cannot be NA.")

  NA.Neg=array(c(S,t,r,dividend,D,k,df))
  NA.Neg.Str=c("S","t","r","dividend","D","k","df")
  app=apply(NA.Neg,1,is.na)
  #Positive
  na.s=which(app==F & NA.Neg<=0)
  if(length(na.s)>0) {errs=paste(NA.Neg.Str[na.s],collapse=" & ")
  stop(cat("Error: '",errs, "' must be positive real number(s).\n"))}
  #Infinite
  na.s2=which(app==F & NA.Neg==Inf)
  if(length(na.s2)>0) {errs=paste(NA.Neg.Str[na.s2],collapse=" & ")
  stop(cat("Error: '",errs, "' cannot be infinite.\n"))}
  #Logical
  stopifnot(is.logical(plot))

  if(df != round(df)) stop("df must be an integer.")
  if(position != "short" & position != "long") stop("Invalid position.")
  if(div.structure != "none" & div.structure != "continuous" & div.structure != "discrete") stop("Invalid dividend structure.")
  if(div.structure=="continuous" & is.na(D)) stop("Missing D.")
  if(div.structure=="continuous" & is.na(D)==FALSE & D <=0) stop("D <=0")
  if((div.structure=="none" | div.structure=="continuous") & is.na(dividend)==FALSE ) warning("Dividend value will be ignored.")
  if((div.structure=="none" | div.structure=="discrete") & is.na(D)==FALSE ) warning("D value will be ignored.")
  if((div.structure=="none" | div.structure=="continuous") & is.na(k)==FALSE ) warning("k value will be ignored.")

  if(div.structure=="discrete"){
    if(is.na(dividend)) stop("Dividend not specified.")
    if(is.na(dividend)==F & dividend<0) stop("Dividend < 0")
    if(is.na(dividend)==F & dividend==0) stop("Set div.structure = none")
    eff.i=exp(r)-1
    int=(1+eff.i)^(1/df)-1
    t1=t*df
    if(!is.na(k)){
      if(int==k){
        price=S*exp(r*t)-dividend*t1/(1+int)*exp(r*t)
      }
      if(int != k){
        price=S*exp(r*t)-dividend*(1-((1+k)/(1+int))^t1)/(int-k)*exp(r*t)
      }
    }
    if(is.na(k)){
      price=S*exp(r*t)-dividend*((1+int)^t1-1)/int
    }
  }

  if(div.structure=="none") price=S*exp(r*t)
  if(div.structure=="continuous") price=S*exp((r-D)*t)

  stock=seq(0,price,length.out=6)
  stock=c(stock,round(seq(price,price*2,length.out=6)))
  stock=sort(unique(round(stock)))
  if(position=="long") payoff=stock-price else payoff=price-stock
  if(plot==T){
    if(position=="long") m="Long Forward\nPayoff" else m="Short Forward\nPayoff"
    plot(stock,payoff,type="l",xlab="Stock Price",main=m,ylab="Payoff",
         ylim=c(min(payoff),max(payoff)),xaxt='n',yaxt='n',col="steelblue",lwd=2)
    abline(h=0,lty=2,col="gray")
    x=round(seq(min(payoff),max(payoff),length.out=8))
    axis(2,at=x,labels=x,las=2)
    axis(1,at=stock) }

  out1=data.frame(stock,payoff)
  names(out1)=c("Stock Price","Payoff")
  out=list(Payoff=out1,Price=price)
  return(out)
}
