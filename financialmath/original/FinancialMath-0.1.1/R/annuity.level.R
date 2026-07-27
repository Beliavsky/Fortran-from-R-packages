annuity.level=function(pv=NA,fv=NA,n=NA,pmt=NA,i=NA,ic=1,pf=1,imm=TRUE,plot=FALSE){

  all=list(pv,fv,n,pmt,i,ic,pf,imm,plot)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(all,length) != 1)==T) stop("All inputs must be of length 1.")
  #Numeric
  num2=list(pv,fv,i,pmt,ic,n,pf)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F)) stop("pv, fv, i, pmt, ic, n, and pf must be numeric.")
  #NA
  nalist=list(ic,pf,imm,plot)
  if(any(lapply(nalist,is.na)==T)) stop("ic, pf, imm, and plot cannot be NA.")

  num=c(pv,fv,i,pmt,ic,pf,n)
  NA.Neg=array(num)
  NA.Neg.Str=c("pv","fv","i","pmt","ic","pf","n")
  app=apply(NA.Neg,1,is.na)
  #Positive
  na.s=which(app==F & NA.Neg<=0)
  if(length(na.s)>0) {
    errs=paste(NA.Neg.Str[na.s],collapse=" & ")
    stop(cat("Error: '",errs, "' must be positive real number(s).\n"))}
  #Infinite
  na.s2=which(app==F & NA.Neg==Inf)
  if(length(na.s2)>0) {errs=paste(NA.Neg.Str[na.s2],collapse=" & ")
  stop(cat("Error: '",errs, "' cannot be infinite.\n"))}
  #Logical
  stopifnot(is.logical(imm),is.logical(plot))

  test=c(n,pmt,i)
  test2=c(pv,fv)
  t1=length(test[which(is.na(test))])
  t2=length(test2[which(is.na(test2))])
  if(t2 == 0) stop("Cannot specify both pv and fv.")
  if(t1 >= 2) stop("Too many unknowns.")
  if(t1==0 & t2==0) stop("No unknowns specified.")
  if(t1==1 & t2==2) stop("Too many unknowns.")
  if(t1==0 & t2==1) stop("One of n, pmt, or i must be NA (unknown) if either pv or fv is known.")

  if(!is.na(n) & n !=round(n)) stop("n must be a positive integer.")

  nom1=NA;nom2=NA
  if(is.na(n)) sn=T else sn=F

  if(is.na(i)){
    if(imm==T) {
      if(is.na(pv)) {i=round(Re(polyroot(c(-fv/pmt+1,rep(1,n-1))))-1,8)
      i=i[which(round(pmt*((1+i)^n-1)/i,2)==round(fv,2))]}
      if(is.na(fv)) {i=round(1/Re(polyroot(c(-pv/pmt,rep(1,n))))-1,8)
      i=i[which(round(pmt*(1-(1+i)^-n)/i,2)==round(pv,2))]}
    }
    if(imm==F) {
      if(is.na(pv)) {i=round(Re(polyroot(c(-fv/pmt,rep(1,n))))-1,8)
      i=i[which(round(pmt*(1+i)*((1+i)^n-1)/i,2)==round(fv,2))]}
      if(is.na(fv)) {i=round(1/Re(polyroot(c(-pv/pmt+1,rep(1,n-1))))-1,8)
      i=i[which(round(pmt*(1+i)*(1-(1+i)^-n)/i,2)==round(pv,2))]}
    }
    if(length(i)>1) {if(any(i>0)) i=i[which(i>0)] else i=i[1];i=i[1]}
    eff.i=(1+i)^pf-1
    n.i=(1+eff.i)^(1/ic)-1
    if(ic != 1) nom1=n.i*ic
    if(pf != 1 & pf != ic) nom2=i*pf
    i=n.i*ic
  }

  eff.i=(1+i/ic)^(ic)-1
  int=(1+eff.i)^(1/pf)-1
  if(imm==T) r=1 else r=1+int
  if(ic !=1) nom1=((1+eff.i)^(1/ic)-1)*ic
  if(pf != ic & pf !=1) nom2=((1+eff.i)^(1/pf)-1)*pf

  if(is.na(pmt)){
    if(is.na(fv)) pmt=pv/((1-(1+int)^-(n))/(int)*r)
    if(is.na(pv)) pmt=fv/(((1+int)^(n)-1)/(int)*r)
  }

  if(is.na(n)){
    if(is.na(fv)) {
    if(round((1-pv/pmt*int/r),2)<0) stop("n is negative.")
    if(round((1-pv/pmt*int/r),2)>0) n=log(1-pv/pmt*int/r)/log(1/(1+int))
    if(round((1-pv/pmt*int/r),2)==0) n=Inf}
    if(is.na(pv)) n=log(fv/pmt*int/r+1)/log(1+int)
  }

  if(is.na(pv)) pv=r*pmt*(1-(1+int)^(-n))/int
  if(is.na(fv)) fv=r*pmt*((1+int)^(n)-1)/int

  if(plot==T){
    if(sn==F){
      plot.new()
      plot(0,0,type="n",axes=F,ann=F,xlim = c(0,ceiling(n)+1),ylim=c(0, 10))
      axis(1,at=seq(0,ceiling(n)),labels=seq(0,ceiling(n)),line=-8)
      if(imm==T){
        text(seq(1,ceiling(n)),rep(5.5,ceiling(n)),labels=rep(round(pmt,2),ceiling(n)),cex=.75) }
      if(imm==F){
        text(seq(0,ceiling(n)-1),rep(5.5,ceiling(n)),labels=rep(round(pmt,2),ceiling(n)),cex=.75) }
      if(pf != 1) axis(1,at=seq(0,ceiling(n),by=pf),labels=seq(0,ceiling(n),by=pf)/pf,line=-5,col="blue")
      text(.05*n,7.2,labels=round(pv,2),cex=.85)
      text(.05*n,7.7,labels="PV")
      text(n-.05*n,7.2,labels=round(fv,2),cex=.85)
      text(n-.05*n,7.7,labels="FV")
      text(n/2,10,labels="Time Diagram",cex=1.2)
      text(n/2,9.4,labels="Level Annuity",cex=1)
      text(n*.4,8.5,bquote("Eff Rate "== .(round(eff.i,4))),cex=.85)
      if(pf != 1) text(n*.4,7.8,bquote( i^(.(pf))== .(round(int*pf,4))),cex=.85)
      i.ic=((1+eff.i)^(1/ic)-1)*ic
      if(ic != 1 & ic != pf) text(n*.4,7.1,bquote( i^(.(ic))== .(round(i.ic,4))),cex=.85)
      if(pf==1) legend(0,1,legend="Years",lty=1,bty="n")
      else legend(0,1,legend=c("Periods","Years"),lty=1,col=c("black","blue"),bty="n",ncol=2)
    } else warning("No time diagram is provided when solving for n.\n")  }

  if(pf==1) {years=n;n=NA} else years=n/pf
  out=c(pv,fv,pmt,eff.i,nom1,nom2,n,years)
  m.out=matrix(out,nrow=length(out))
  rownames(m.out)=c("PV","FV","PMT","Eff Rate",paste("i^(",round(ic,2),")",sep=""),
                    paste("i^(",round(pf,2),")",sep=""),"Periods","Years")
  na=apply(m.out, 1, function(x) all(is.na(x)))
  m.out=as.matrix(m.out[!na,])
  colnames(m.out)= if(is.infinite(n) | is.infinite(years)) "Level Perpetuity" else "Level Annuity"
  return(m.out)
}
