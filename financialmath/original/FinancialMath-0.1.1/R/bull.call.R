bull.call=function(S,K1,K2,r,t,price1,price2,plot=FALSE){

  all=list(S,K1,K2,r,t,price1,price2,plot)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(all,length) != 1)==T) stop("All inputs must be of length 1.")
  #Numeric
  num2=list(S,K1,K2,r,t,price1,price2)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F)) stop("S K1, K2, r, t, price1, and price2 must be numeric.")
  #NA
  if(any(lapply(all,is.na)==T)) stop("Cannot input any variables as NA.")

  NA.Neg=c(r,S,K1,K2,t,price1,price2)
  #Positive
  if(any(NA.Neg<=0)) stop("All numeric variables must be positive.")
  #Infinite
  if(any(NA.Neg==Inf)) stop("Cannot input any variables as infinite.")
  #Logical
  stopifnot(is.logical(plot))

  if(K1>=S) stop("K1 must be less than S.")
  if(S>=K2) stop("K2 must be greater than S.")

  #long
  callP1 = price1
  #short
  callP2 = price2

  stock=unique(round(seq(0,K1,length.out=6)))
  stock=c(stock,round(seq(K1,K2,length.out=4)))
  stock=c(stock,round(seq(K2,K2+K1,length.out=6)))
  stock=unique(stock)
  stock2=unique(round(c(seq(0,K1,length.out=6),seq(K2,K2+K1,length.out=6))))
  payoff=rep(0,length(stock))
  profit=rep(0,length(stock))

  for(i in 1:length(stock)){
    if(stock[i]<=K1) payoff[i]=0
    if(stock[i]>=K2) payoff[i]=K2-K1
    if(stock[i]<K2 & stock[i]>K1) payoff[i]=stock[i]-K1
    profit[i]=payoff[i]+(callP2-callP1)*exp(r*t)
  }

  if(plot==T){
    plot(stock,profit,type="l",xlab="Stock Price",main="Bull Call Spread\nPayoff and Profit",ylab="$",
         ylim=c(min(profit,payoff),max(profit,payoff)),xaxt='n',yaxt='n',col="steelblue",lwd=2)
    lines(stock,payoff,lty=2,lwd=2,col="firebrick")
    abline(h=0,lty=2,col="gray")
    y=round(seq(min(payoff,profit),max(payoff,profit),length.out=8))
    axis(2,at=y,labels=y,las=2)
    axis(1,at=if(K2-K1<7) stock2[-which(stock2==K2)] else stock2,las=2)
    legend("topleft",c("Profit","Payoff"),lty=c(1,2),col=c("steelblue","firebrick"),lwd=c(2,2)) }

  out1=data.frame(stock,payoff,profit)
  names(out1)=c("Stock Price","Payoff","Profit")
  out2=matrix(c(callP1,callP2,callP1-callP2),nrow=3)
  rownames(out2)=c("Call (K1)","Call (K2)","Net Cost")
  colnames(out2)=c("Premiums")
  out=list(Payoff=out1,Premiums=out2)
  return(out)
}
