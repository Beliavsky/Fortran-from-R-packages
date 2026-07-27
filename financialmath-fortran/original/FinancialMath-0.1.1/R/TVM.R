TVM=function(pv=NA,fv=NA,n=NA,i=NA,ic=1,plot=FALSE) {

  all=list(pv,fv,n,i,ic,plot)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(all,length) != 1)==T) stop("All inputs must be of length 1.")
  #Numeric
  num2=list(pv,fv,n,i,ic)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F)) stop("pv, fv, n, i, and ic must be numeric.")
  #NA
  nalist=list(plot,ic)
  if(any(lapply(nalist,is.na)==T)) stop("Cannot input ic or plot as NA.")

  NA.Neg=array(c(pv,fv,n,i,ic))
  NA.Neg.Str=c("pv","fv","n","i","ic")
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

  test=c(pv,fv,n,i)
  if(length(test[which(is.na(test))])>=2) stop("Too many unknowns.")
  if(length(test[which(is.na(test))])==0) stop("Must have one unknown.")
  if(!is.na(fv) & pv>fv & !is.na(pv)) stop("PV cannot be greater than FV.")

  nom1=NA
  if(!is.na(i)) {int=(1+i/ic)^ic-1
    if(ic!=1) nom1=i}

  if(is.na(pv)){pv=fv/(1+int)^(n);solve="PV"}
  if(is.na(fv)){fv=pv*(1+int)^(n);solve="FV"}
  if(is.na(n)){n=log(fv/pv)/log(1+int);solve="n"}
  if(is.na(i)){int=(fv/pv)^(1/n)-1;solve="i"
    if(ic!=1) nom1=((1+int)^(1/ic)-1)*ic}

  if (plot==TRUE) {
    x=unique(sort(c(seq(0,max(n),by=1),n)))
    y=pv*(1+int)^x
    plot.new()
    plot(x,y,ylab="Value",xlab="Period",xlim=c(0,max(x)),
         ylim=c(pv,max(y)),main="Time Value",pch=20)
    lines(x,y,lty=3)
    if(solve=="PV") points(0,pv,pch=20,col="blue")
    if(solve=="FV") points(n,fv,pch=20,col="blue")
  }

  out=c(pv,fv,n,int,nom1)
  m.out=matrix(out,nrow=length(out))
  rownames(m.out)=c("PV","FV","Periods",
                    "Eff Rate",paste("i^(",round(ic,2),")",sep=""))
  na=apply(m.out, 1, function(x) all(is.na(x)))
  m.out=as.matrix(m.out[!na,])
  colnames(m.out)=c("TVM")
  return(m.out)
}

