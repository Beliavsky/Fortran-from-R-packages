annuity.geo=function(pv=NA,fv=NA,n=NA,p=NA,k=NA,i=NA,ic=1,pf=1,imm=TRUE,plot=FALSE){

  all=list(pv,fv,n,p,k,i,ic,pf,imm,plot)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(all,length) != 1)==T) stop("All inputs must be of length 1.")
  #Numeric
  num2=list(pv,fv,n,p,k,i,ic,pf)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F)) stop("pv, fv, n, p, k, i, ic, and pf must be numeric.")
  #NA
  nalist=list(ic,pf,imm,plot)
  if(any(lapply(nalist,is.na)==T)) stop("ic, pf, imm, and plot cannot be NA.")

  NA.Neg=array(c(pv,fv,n,p,i,ic,pf))
  NA.Neg.Str=c("pv","fv","n","p","i","ic","pf")
  app=apply(NA.Neg,1,is.na)
  #Positive
  na.s=which(app==F & NA.Neg<=0)
  if(length(na.s)>0) {errs=paste(NA.Neg.Str[na.s],collapse=" & ")
  stop(cat("Error: '",errs, "' must be positive real number(s).\n"))}
  #Infinite
  na.s2=which(app==F & NA.Neg==Inf)
  if(length(na.s2)>0) {errs=paste(NA.Neg.Str[na.s2],collapse=" & ")
  stop(cat("Error: '",errs, "' cannot be infinite.\n"))}
  if(!is.na(k) & k==Inf) stop("k cannot be infinite.")
  #Logical
  stopifnot(is.logical(imm),is.logical(plot))

  if(!is.na(n) & n !=round(n)) stop("n must be a positive integer.")
  if(is.na(n)) sn=T else sn=F
  nom1=NA;nom2=NA

  test=c(n,i);t1=length(test[which(is.na(test))])
  test2=c(pv,fv);t2=length(test2[which(is.na(test2))])
  test3=c(p,k);t3=length(test3[which(is.na(test3))])
  test4=c(n,i,p,q);t4=length(test4[which(is.na(test4))])
  if(t2==0) stop("Cannot specify both pv and fv.")
  if(t1>=2) stop("Too many unknowns.")
  if(t3>=2) stop("Cannot have both p and k unknown")
  if(t4==0 & t2==0) stop("No unknowns specified.")
  if(t4==1 & t2==2) stop("Too many unknowns.")
  if(t1==0 & t3==0 & t2==1) stop("one of n, p, k, or i must be unknown if either pv or fv is known.")

  if(is.na(k)){
    int=((1+i/ic)^ic)^(1/pf)-1;r=1+int
    if(imm==T){
      if(is.na(pv)){k=round(Re(polyroot(c(-fv+p*(1+int)^(n-1),rev(p*(1+int)^seq(0,n-2)))))-1,8)
      if(all(round(int,3)!=round(k,3))) k=k[which(round(fv,2)==round(p*(1+int)^n*(1-((1+k)/(1+int))^n)/(int-k),2))]
      else k=k[which(round(fv,2)==round(p/(1+int)*n*(1+int)^n,2))]}
      if(is.na(fv)){k=round(Re(polyroot(c(-pv+p/(1+int),p*(1+int)^(-1*seq(2,n)))))-1,8)
      if(all(round(int,3)!=round(k,3))) k=k[which(round(pv,2)==round(p*(1-((1+k)/(1+int))^n)/(int-k),2))]
      else k=k[which(round(pv,2)==round(p/(1+int)*n,2))]}
    }
    if(imm==F){
      if(is.na(pv)){k=round(Re(polyroot(c(-fv+p*(1+int)^n,rev(p*(1+int)^seq(1,n-1)))))-1,8)
      if(all(round(int,3)!=round(k,3))) k=k[which(round(fv,2)==round(p*r*(1+int)^n*(1-((1+k)/(1+int))^n)/(int-k),2))]
      else k=k[which(round(fv,2)==round(p/(1+int)*r*n*(1+int)^n,2))]}
      if(is.na(fv)){k=round(Re(polyroot(c(-pv+p,p*(1+int)^(-1*seq(1,n-1)))))-1,8)
      if(all(round(int,3)!=round(k,3))) k=k[which(round(pv,2)==round(r*p*(1-((1+k)/(1+int))^n)/(int-k),2))]
      else k=k[which(round(pv,2)==round(p/(1+int)*n*r,2))]}
    }
    k=k[1]
  }

  if(is.na(i)){
    if(imm==T) {
      if(is.na(pv)) {i=round(Re(polyroot(c(-fv+p*(1+k)^(n-1),rev(p*(1+k)^seq(0,n-2)))))-1,8)
      if(all(round(i,2)!=round(k,2))) i=i[which(round((p*(1-((1+k)/(1+i))^n)/(i-k))*(1+i)^n,1)==round(fv,1))]}
      if(is.na(fv)) {i=round(1/Re(polyroot(c(-pv,p*(1+k)^seq(0,n-1))))-1,8)
      if(all(round(i,2)!=round(k,2))) i=i[which(round(p*(1-((1+k)/(1+i))^n)/(i-k),1)==round(pv,1))]}
    }
    if(imm==F) {
      if(is.na(pv)) {i=round(Re(polyroot(c(-fv,rev(p*(1+k)^seq(0,n-1)))))-1,8)
      if(all(round(i,2)!=round(k,2))) i=i[which(round(p*(1-((1+k)/(1+i))^(n))/(i-k)*(1+i)^(n+1),1)==round(fv,1))]}
      if(is.na(fv)) {i=round(1/Re(polyroot(c(-pv+p,(p*(1+k)^seq(1,n-1)))))-1,8)
      if(all(round(i,2)!=round(k,2))) i=i[which(round(p*(1-((1+k)/(1+i))^(n))/(i-k)*(1+i),1)==round(pv,1))]}
    }

    if(any(i>0 & i<.1)) i=i[which(i>0 & i<.1)]
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

  if(is.na(p)){
    if(is.na(fv)) {
      if(round(int,3)!=round(k,3)) p=(pv/(r*(1-((1+k)/(1+int))^n)/(int-k)))
      if(round(int,3)==round(k,3)) p=(pv/r)*(1+int)/n
    }
    if(is.na(pv)) {
      if(round(int,3)!=round(k,3)) p=(fv/(r*(1-((1+k)/(1+int))^n)/(int-k))/(1+int)^n)
      if(round(int,3)==round(k,3)) p=(fv/r)*((1+int)/n)*(1+int)^(-n)
    }
  }

  if(is.na(n)){
    if(is.na(fv)) {
      if(round(int,3)!=round(k,3)) n=log(-(pv/p/r*(int-k)-1))/log((1+k)/(1+int))
      if(round(int,3)==round(k,3)) n=(pv/p)/r*(1+int)
    }
    if(is.na(pv)) {
      if(round(int,3)!=round(k,3)) {n=seq(0,1000,by=.001)
      n=n[which(round(fv,2)==round(p*r*(1+int)^n*(1-((1+k)/(1+int))^n)/(int-k),2))]
      }
      if(round(int,3)==round(k,3)) {n=seq(0,1000,by=.001)
      n=n[which(round(n*(1+int)^n,2)==round((fv/p)*(1+int)/r,2))]
      }
      n=median(n)}
  }

  if(is.na(pv)) {
    if(round(int,3)!=round(k,3)) pv=(r*p*(1-((1+k)/(1+int))^n)/(int-k))
    if(round(int,3)==round(k,3)) pv=p/(1+int)*n*r
  }
  if(is.na(fv)) {
    if(round(int,3)!=round(k,3)) fv=(r*p*(1-((1+k)/(1+int))^n)/(int-k))*(1+int)^n
    if(round(int,3)==round(k,3)) fv=p/(1+int)*n*r*(1+int)^n
  }

  if(plot==T){
    if(sn==F){
      plot.new()
      plot(0,0,type="n",axes=F,ann=F,xlim = c(0,ceiling(n)+1),ylim=c(0, 10))
      axis(1,at=seq(0,ceiling(n)),labels=seq(0,ceiling(n)),line=-8)
      if(imm==T){
        text(seq(1,ceiling(n)),rep(5.5,ceiling(n)),labels=round(p*(1+k)^seq(0,n-1),2),cex=.7)
      }
      if(imm==F){
        text(seq(0,ceiling(n)-1),rep(5.5,ceiling(n)),labels=round(p*(1+k)^seq(0,n-1),2),cex=.7)
      }
      if(pf != 1) axis(1,at=seq(0,ceiling(n),by=pf),labels=seq(0,ceiling(n),by=pf)/pf,line=-5,col="blue")
      text(.05*n,7.2,labels=round(pv,2),cex=.85)
      text(.05*n,7.7,labels="PV")
      text(n-.05*n,7.2,labels=round(fv,2),cex=.85)
      text(n-.05*n,7.7,labels="FV")
      text(n/2,10,labels="Time Diagram",cex=1.2)
      if(round(k,2)==0) m="Level Annuity" else m="Geometric Annuity"
      text(n/2,9.4,labels=m,cex=1)
      text(n*.4,8.5,bquote("Eff Rate "== .(round(eff.i,4))),cex=.85)
      if(pf != 1) text(n*.4,7.8,bquote( i^(.(pf))== .(round(int*pf,4))),cex=.85)
      i.ic=((1+eff.i)^(1/ic)-1)*ic
      if(ic != 1 & ic != pf) text(n*.4,7.1,bquote( i^(.(ic))== .(round(i.ic,4))),cex=.85)
      if(pf==1) legend(0,1,legend="Years",lty=1,bty="n") else legend(0,1,legend=c("Periods","Years"),lty=1,col=c("black","blue"),bty="n",ncol=2)
    } else warning("No time diagram is provided when solving for n.\n")  }

  if(pf==1) {years=n;n=NA} else years=n/pf
  out=c(pv,fv,p,k,eff.i,nom1,nom2,n,years)
  m.out=matrix(out,nrow=length(out))
  rownames(m.out)=c("PV","FV","P","K","Eff Rate",paste("i^(",round(ic,2),")",sep=""),
                    paste("i^(",round(pf,2),")",sep=""),"Periods","Years")
  na=apply(m.out, 1, function(x) all(is.na(x)))
  m.out=as.matrix(m.out[!na,])
  if(round(k,2)==0) m2="Level Annuity" else m2="Geometric Annuity"
  colnames(m.out)=m2
  return(m.out)
}
