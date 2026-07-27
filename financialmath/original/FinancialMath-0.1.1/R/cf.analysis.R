cf.analysis=function(cf,times,i,plot=FALSE,time.d=FALSE){

  all=list(i,cf,times,plot,time.d)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(list(i,plot,time.d),length) != 1)==T) stop("i, plot, and time.d must be of length 1.")
  #Numeric
  if(!is.vector(cf) | !is.numeric(cf)) stop("cf must be a numeric vector.")
  if(!is.vector(times) | !is.numeric(times)) stop("times must be a numeric vector.")
  num2=list(i)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F)) stop("i must be numeric.")
  #NA
  if(any(is.na(cf)) | any(is.na(times)) | any(is.na(c(i,plot,time.d)))) stop("Cannot input any variables as NA.")
  #Infinite
  if(i==Inf) stop("i cannot be infinite.")
  if(any(cf==Inf)) stop("Cannot have infinite in cf.")
  if(any(times==Inf)) stop("Cannot have infinite in times.")
  #Positive
  if(i<0) stop("i must be positive.")
  if(any(times<0)) stop("Cannot have negative values in times.")
  #Logical
  stopifnot(is.logical(plot),is.logical(time.d))


  if(length(cf)==0 | length(times)==0 ) stop("Missing cf or times information.")
  if(length(cf) != length(times)) stop("Number of payments not equal to number of period values.")

  d=unique(times[which(duplicated(times)==T)])
  if(length(d)>0){ cfd=rep(0,length(d))
  for(r in 1:length(d)) cfd[r]=sum(cf[which(times==d[r])])}

  p=rep(0,max(times))
  p[times]=cf
  if(length(d)>0) {for(r in 1:length(d)) p[d[r]]=cfd[r]}
  pv=0
  for(m in 1:length(times)) pv=pv+cf[m]/(1+i)^times[m]

  macd.n=0
  for(m in 1:length(times)) macd.n=macd.n+(times[m]*cf[m]/(1+i)^times[m])
  if(pv==0) macd=0 else macd=macd.n/pv
  modd.n=0
  for(m in 1:length(times)) modd.n=modd.n+(times[m]*cf[m]/(1+i)^(times[m]+1))
  if(pv==0) modd=0 else modd=modd.n/pv
  macc.n=0
  for(m in 1:length(times)) macc.n=macc.n+(times[m]^2*cf[m]/(1+i)^times[m])
  if(pv==0) macc=0 else macc=macc.n/pv
  modc.n=0
  for(m in 1:length(times)) modc.n=modc.n+(times[m]*(times[m]+1)*cf[m]/(1+i)^(times[m]+2))
  if(pv==0) modc=0 else modc=modc.n/pv


  if(plot==T){
    int=seq(0,2*i,by=.001)
    pv.g=rep(0,length(int))
    for(x in 1:length(int)){
      for(m in 1:length(times)) pv.g[x]=pv.g[x]+cf[m]/(1+int[x])^times[m] }
    plot(int,pv.g,type="l",ylab="Present Value",xlab="Interest Rate (per period)",main="Convexity")
    abline(h=pv,col="gray",lty=2)
    abline(v=i,col="gray",lty=2)
  }

  if(time.d==T){
    if(plot==T) dev.new() else plot.new()
    x.pv=p[which(p!=0)]
    if(length(x.pv)==0) x.pv=0
    x.t=unique(times)[order(unique(times))]
    plot(0,0,type="n",axes=F,ann=F,xlim = c(-1,max(x.t)+1),ylim=c(0, 10))
    axis(1,at=c(0,x.t),labels=c(0,x.t),line=-5)
    par(mgp=c(3,2,0))
    par(col.axis="blue")
    text(0,6.5,labels=round(pv,2),cex=.8,col="blue")
    par(col.axis="black")
    text(x.t,rep(5.5,max(x.t)),labels=x.pv,cex=.85)
    text((max(x.t))/2,10,labels="Time Diagram",cex=1.2)
    text((max(x.t))/2,8.5,bquote(i== .(round(i,4))),cex=.85)
    legend(0,1,legend="Period",lty=0,bty="n")
    par(mgp=c(3,1,0))
  }

  out=matrix(c(pv,macd,modd,macc,modc),nrow=5)
  rownames(out)=c("PV","MAC D","MOD D","MAC C","MOD C")
  na=apply(out, 1, function(x) all(is.na(x)))
  out=as.matrix(out[!na,])
  colnames(out)="CF Analysis"
  return(out)
}
