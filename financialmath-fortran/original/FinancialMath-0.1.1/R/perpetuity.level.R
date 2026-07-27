perpetuity.level=function(pv=NA,pmt=NA,i=NA,ic=1,pf=1,imm=TRUE){

  all=list(pv,pmt,i,ic,pf,imm)
  #Null
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(all,length) != 1)==T) stop("All inputs must be of length 1.")
  #Numeric
  num2=list(pv,pmt,i,ic,pf)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F)) stop("pv, pmt, i, ic, and pf must be numeric.")
  #NA
  nalist=list(ic,pf,imm)
  if(any(lapply(nalist,is.na)==T)) stop("ic, pf, and imm cannot be NA.")
  #Logical
  stopifnot(is.logical(imm))

  NA.Neg=array(c(pv,pmt,i,ic,pf))
  NA.Neg.Str=c("pv","pmt","i","ic","pf")
  app=apply(NA.Neg,1,is.na)
  #Positive
  na.s=which(app==F & NA.Neg<=0)
  if(length(na.s)>0) {errs=paste(NA.Neg.Str[na.s],collapse=" & ")
  stop(cat("Error: '",errs, "' must be positive real number(s).\n"))}
  #Infinite
  na.s2=which(app==F & NA.Neg==Inf)
  if(length(na.s2)>0) {errs=paste(NA.Neg.Str[na.s2],collapse=" & ")
  stop(cat("Error: '",errs, "' cannot be infinite.\n"))}

  test=c(pv,pmt,i)
  if(length(test[which(is.na(test))]) !=1) stop("Must have exactly 1 NA variable.")
  nom1=NA;nom2=NA

  if(is.na(i)){
    if(imm==T) int=pmt/pv
    if(imm==F) int=(-pmt/pv+1)^(-1)-1
    eff.i=(1+int)^pf-1
    i=((1+eff.i)^(1/ic)-1)*ic
    if(ic != 1) nom1=i
    if(pf != 1 & pf != ic) nom2=int*pf
  }

  eff.i=(1+i/ic)^ic-1
  int=(1+eff.i)^(1/pf)-1
  if(imm==T) r=1 else r=1+int
  if(is.na(pmt)) pmt=pv/r*int
  if(is.na(pv)) pv=pmt*r/int

  if(ic !=1) nom1=((1+eff.i)^(1/ic)-1)*ic
  if(pf != ic & pf !=1) nom2=((1+eff.i)^(1/pf)-1)*pf

  out=matrix(c(pv,pmt,eff.i,nom1,nom2),nrow=5)
  rownames(out)=c("PV","PMT","Eff Rate",paste("i^(",round(ic,2),")",sep=""),
                  paste("i^(",round(pf,2),")",sep=""))
  na=apply(out, 1, function(x) all(is.na(x)))
  out=as.matrix(out[!na,])
  colnames(out)=c("Level Perpetuity")
  return(out)
}
