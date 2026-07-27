strangle = function(S,K1,K2,r,t,price1,price2,plot=FALSE){

  all=list(S,K1,K2,r,t,price1,price2,plot)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(all,length) != 1)==T) stop("All inputs must be of length 1.")
  #Numeric
  num2=list(S,K1,K2,r,t,price1,price2)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F)) stop("S, K1, K2, r, t, price1, and price2 must be numeric.")
  #NA
  nalist=list(S,K1,K2,r,t,price1,price2,plot)
  if(any(lapply(nalist,is.na)==T)) stop("Cannot input any variables as NA.")
  #Logical
  stopifnot(is.logical(plot))

  NA.Neg=array(c(S,K1,K2,r,t,price1,price2))
  NA.Neg.Str=c("S","K1","K2","r","t","price1","price2")
  app=apply(NA.Neg,1,is.na)
  #Positive
  na.s=which(app==F & NA.Neg<=0)
  if(length(na.s)>0) {errs=paste(NA.Neg.Str[na.s],collapse=" & ")
  stop(cat("Error: '",errs, "' must be positive real number(s).\n"))}
  #Infinite
  na.s2=which(app==F & NA.Neg==Inf)
  if(length(na.s2)>0) {errs=paste(NA.Neg.Str[na.s2],collapse=" & ")
  stop(cat("Error: '",errs, "' cannot be infinite.\n"))}

  if(K1>=K2) stop("K1 must be less than K2.")
  if(K2<=S | K1>=S) stop("K2 must be > S and K1 must be < S")

  putP = price1
  callP = price2

  stock=unique(round(seq(0,K1,length.out=6)))
  stock=c(stock,round(seq(K1,K2,length.out=4)))
  stock=c(stock,round(seq(K2,K2+K1,length.out=6)))
  stock=unique(stock)
  stock2=unique(round(c(seq(0,K1,length.out=6),seq(K2,K2+K1,length.out=6))))
  payoff=rep(0,length(stock))
  profit=rep(0,length(stock))
  for(i in 1:length(stock)){
    if(stock[i]<=K1) payoff[i]=K1-stock[i]
    if(stock[i]>=K2) payoff[i]=stock[i]-K2
    if(stock[i]<K2 & stock[i]>K1) payoff[i]=0
    profit[i]=payoff[i]+(-callP-putP)*exp(r*t)
  }

  if(plot==T){
    plot(stock,profit,type="l",xlab="Stock Price",main="Strangle\nPayoff and Profit",ylab="$",
         ylim=c(min(profit,payoff),max(profit,payoff)),xaxt='n',yaxt='n',col="steelblue",lwd=2)
    lines(stock,payoff,lty=2,lwd=2,col="firebrick")
    abline(h=0,lty=2,col="gray")
    y=round(seq(min(payoff,profit),max(payoff,profit),length.out=8))
    axis(2,at=y,labels=y,las=2)
    axis(1,at=if(K2-K1<7) stock2[-which(stock2==round(K2))] else stock2,las=2)
    legend("top",c("Profit","Payoff"),lty=c(1,2),col=c("steelblue","firebrick"),lwd=c(2,2)) }


  out1=data.frame(stock,payoff,profit)
  names(out1)=c("Stock Price","Payoff","Profit")
  out2=matrix(c(putP,callP,callP+putP),nrow=3)
  rownames(out2)=c("Put (K1)","Call (K2)","Net Cost")
  colnames(out2)=c("Premiums")
  out=list(Payoff=out1,Premiums=out2)
  return(out)
}
