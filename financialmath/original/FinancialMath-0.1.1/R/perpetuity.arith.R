perpetuity.arith=function(pv=NA,p=NA,q=NA,i=NA,ic=1,pf=1,imm=TRUE){

  all=list(pv,p,q,i,ic,pf,imm)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(all,length) != 1)==T) stop("All inputs must be of length 1.")
  #Numeric
  num2=list(pv,p,i,ic,pf,q)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F)) stop("pv, p, q, i, ic, and pf must be numeric.")
  #NA
  nalist=list(ic,pf,imm)
  if(any(lapply(nalist,is.na)==T)) stop("Cannot input ic, pf, and imm as NA.")
  #Logical
  stopifnot(is.logical(imm))

  num=c(pv,p,i,ic,pf)
  NA.Neg=array(num)
  NA.Neg.Str=c("pv","p","i","ic","pf")
  app=apply(NA.Neg,1,is.na)
  #Positive
  na.s=which(app==F & NA.Neg<=0)
  if(length(na.s)>0) {errs=paste(NA.Neg.Str[na.s],collapse=" & ")
  stop(cat("Error: '",errs, "' must be positive real number(s).\n"))}
  if(q<0) stop("q cannot be negative.")
  #Infinite
  na.s2=which(app==F & NA.Neg==Inf)
  if(length(na.s2)>0) {errs=paste(NA.Neg.Str[na.s2],collapse=" & ")
  stop(cat("Error: '",errs, "' cannot be infinite.\n"))}

  test=c(pv,p,q,i)
  if(length(test[which(is.na(test))]) !=1) stop("Must have exactly 1 NA/unknown variable.")
  nom1=NA;nom2=NA

  if(is.na(i)){
    if(imm==T) i=1/Re(polyroot(c(-pv,p,q)))
    if(imm==F) i=1/Re(polyroot(c(-pv+p,p+q,q)))
    if(length(i)>1) {i=i[which(i>0)];i=i[1]}
    eff.i=(1+i)^pf-1
    n.i=(1+eff.i)^(1/ic)-1
    if(ic != 1) nom1=n.i*ic
    if(pf != 1 & pf != ic) nom2=i*pf
    i=n.i*ic
  }

  eff.i=(1+i/ic)^ic-1
  int=(1+eff.i)^(1/pf)-1
  if(imm==T) r=1 else r=1+int
  if(ic !=1) nom1=((1+eff.i)^(1/ic)-1)*ic
  if(pf != ic & pf !=1) nom2=((1+eff.i)^(1/pf)-1)*pf

  if(is.na(p)) p=((pv/r-q/(int^2))*int)
  if(is.na(q)) q=(pv/r-p/int)*int^2
  if(is.na(pv)) pv=(p/int+q/(int^2))*r

  out=c(pv,p,q,eff.i,nom1,nom2)
  m.out=matrix(out,nrow=length(out))
  rownames(m.out)=c("PV","P","Q","Eff Rate",paste("i^(",round(ic,2),")",sep=""),
                    paste("i^(",round(pf,2),")",sep=""))
  na=apply(m.out, 1, function(x) all(is.na(x)))
  m.out=as.matrix(m.out[!na,])
  if(round(q,2)==0) m2="Level Perpetuity" else m2="Arithmetic Perpetuity"
  colnames(m.out)=m2
  return(m.out)
}
