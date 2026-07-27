protective.put=function(S,K,r,t,sd,price=NA,plot=FALSE){

  all=list(S,K,r,t,sd,plot,price)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(all,length) != 1)==T) stop("All inputs must be of length 1.")
  #Numeric
  num2=list(S,K,r,t,sd,price)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F)) stop("S, K, r, t, sd, and price must be numeric.")
  #NA
  nalist=list(S,K,r,t,sd,plot)
  if(any(lapply(nalist,is.na)==T)) stop("Cannot input any variables, but price, as NA.")
  #Logical
  stopifnot(is.logical(plot))

  NA.Neg=array(c(S,K,r,t,sd,price))
  NA.Neg.Str=c("S","K","r","t","sd","price")
  app=apply(NA.Neg,1,is.na)
  #Positive
  na.s=which(app==F & NA.Neg<=0)
  if(length(na.s)>0) {errs=paste(NA.Neg.Str[na.s],collapse=" & ")
  stop(cat("Error: '",errs, "' must be positive real number(s).\n"))}
  #Infinite
  na.s2=which(app==F & NA.Neg==Inf)
  if(length(na.s2)>0) {errs=paste(NA.Neg.Str[na.s2],collapse=" & ")
  stop(cat("Error: '",errs, "' cannot be infinite.\n"))}

  if(is.na(price)){
    d1 = (log(S/K)+(r+sd^2/2)*t)/(sd*sqrt(t))
    d2 = d1 - sd * sqrt(t)
    putP = K*exp(-r*t) * pnorm(-d2) - S*pnorm(-d1)} else putP=price

  stock=seq(0,K,length.out=6)
  stock=c(stock,round(seq(K,K*2,length.out=6)))
  stock=unique(round(stock))
  payoff=rep(0,length(stock))
  profit=rep(0,length(stock))

  for(i in 1:length(stock)){
    if(stock[i]<=K) payoff[i]=K-S
    if(stock[i]>K) payoff[i]=stock[i]-S
    profit[i]=payoff[i]-(putP)*exp(r*t)
  }

  if(plot==T){
    lpos="topleft";m="Pretective Put\nPayoff and Profit"
    plot(stock,profit,type="l",xlab="Stock Price",main=m,ylab="$",
         ylim=c(min(profit,payoff),max(profit,payoff)),xaxt='n',yaxt='n',lwd=2,col="steelblue")
    lines(stock,payoff,lty=2,lwd=2,col="firebrick")
    abline(h=0,lty=2,col="gray")
    y=round(seq(min(payoff,profit),max(payoff,profit),length.out=8))
    axis(2,at=y,labels=y,las=2)
    axis(1,at=stock)
    legend(lpos,c("Profit","Payoff"),lty=c(1,2),col=c("steelblue","firebrick"),lwd=c(2,2)) }

  out1=data.frame(stock,payoff,profit)
  names(out1)=c("Stock Price","Payoff","Profit")
  out=list(Payoff=out1,Premium=putP)
  return(out)
}
