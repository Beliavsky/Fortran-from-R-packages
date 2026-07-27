straddle=function(S,K,r,t,price1,price2,position,plot=FALSE){

  all=list(S,K,r,t,price1,price2,plot,position)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(all,length) != 1)==T) stop("All inputs must be of length 1.")
  #Numeric
  num2=list(S,K,r,t,price1,price2)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F)) stop("S, K, r, t, price1, and price2 must be numeric.")
  #NA
  nalist=list(S,K,r,t,price1,price2,plot,position)
  if(any(lapply(nalist,is.na)==T)) stop("Cannot input any variables, as NA.")
  #Logical
  stopifnot(is.logical(plot))

  NA.Neg=array(c(S,K,r,t,price1,price2))
  NA.Neg.Str=c("S","K","r","t","price1","price2")
  app=apply(NA.Neg,1,is.na)
  #Positive
  na.s=which(app==F & NA.Neg<=0)
  if(length(na.s)>0) {errs=paste(NA.Neg.Str[na.s],collapse=" & ")
  stop(cat("Error: '",errs, "' must be positive real number(s).\n"))}
  #Infinite
  na.s2=which(app==F & NA.Neg==Inf)
  if(length(na.s2)>0) {errs=paste(NA.Neg.Str[na.s2],collapse=" & ")
  stop(cat("Error: '",errs, "' cannot be infinite.\n"))}

  if(position != "long" & position != "short") stop("Position must be either short or long.")

  callP = price1
  putP = price2

  stock=seq(0,K,length.out=6)
  stock=c(stock,round(seq(K,K*2,length.out=6)))
  stock=unique(round(stock))
  payoff=rep(0,length(stock))
  profit=rep(0,length(stock))
  if (position=="long"){
    for(i in 1:length(stock)){
      if(stock[i]==K) payoff[i]=0
      if(stock[i]<K) payoff[i]=K-stock[i]
      if(stock[i]>K) payoff[i]=stock[i]-K
      profit[i]=payoff[i]-(callP+putP)*exp(r*t)
    }
    cost=callP+putP
  }
  if (position=="short"){
    for(i in 1:length(stock)){
      if(stock[i]==K) payoff[i]=0
      if(stock[i]<K) payoff[i]=stock[i]-K
      if(stock[i]>K) payoff[i]=K-stock[i]
      profit[i]=payoff[i]+(callP+putP)*exp(r*t)
    }
    cost=-(callP+putP)
  }

  if(plot==T){
    if(position=="long") {lpos="bottomright";m="Purchased Straddle\nPayoff and Profit"} else {lpos="topright";m="Written Straddle\nPayoff and Profit"}
    plot(stock,profit,type="l",xlab="Stock Price",main=m,ylab="$",
         ylim=c(min(profit,payoff),max(profit,payoff)),xaxt='n',yaxt='n',col="steelblue",lwd=2)
    lines(stock,payoff,lty=2,lwd=2,col="firebrick")
    abline(h=0,lty=2,col="gray")
    x=round(seq(min(payoff,profit),max(payoff,profit),length.out=8))
    axis(2,at=x,labels=x,las=2)
    axis(1,at=stock)
    legend(lpos,c("Profit","Payoff"),lty=c(1,2),col=c("steelblue","firebrick"),lwd=c(2,2)) }

  out1=data.frame(stock,payoff,profit)
  names(out1)=c("Stock Price","Payoff","Profit")
  out2=matrix(c(callP,putP,cost),nrow=3)
  rownames(out2)=c("Call","Put","Net Cost")
  colnames(out2)=c("Premiums")
  out=list(Payoff=out1,Premiums=out2)
  return(out)
}
