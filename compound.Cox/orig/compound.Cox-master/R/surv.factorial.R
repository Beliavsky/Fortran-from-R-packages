surv.factorial=function(t.vec,d.vec,group,copula=CG.Clayton,alpha,
                    R=1000,t.upper=min(tapply(t.vec,group,max)),
                    C=NULL,S.plot=TRUE,mark.time=FALSE){

  d=length(levels(factor(group)))
  group=as.factor(group)
  A=diag(rep(1,d))%x%t(rep(1,d))/d

  p_func=function(t.vec,d.vec,group){

    w.mat=matrix(0,d,d)
    for(i in 1:d){
      ni=sum(group==i)
      ti=t.vec[group==i]
      di=d.vec[group==i]
      CGi=copula(ti,di,alpha,S.plot=F)
      Si.pos=CGi$surv
      Si.pos[ni]=0
      Si.upper=Si.pos[sum(CGi$time<=t.upper)]
      Si.neg=c(1,CGi$surv[-ni])
      for(l in 1:d){
        nl=sum(group==l)
        tl=t.vec[group==l]
        dl=d.vec[group==l]
        CGl=copula(tl,dl,alpha,S.plot=F)
        Sl.pos=CGl$surv
        Sl.pos[nl]=0
        Sl.upper=Sl.pos[sum(CGl$time<=t.upper)]
        Sl.neg=c(1,CGl$surv[-nl])
        dSl=Sl.neg-Sl.pos
        dSl[CGl$time>t.upper]=0
        temp=colSums(matrix(sort(ti),ni,nl)<=matrix(sort(tl),ni,nl,byrow=T))
        temp=pmax(temp,1)
        w.mat[i,l]=sum((Si.pos[temp]+Si.neg[temp])*dSl/2)+Si.upper*Sl.upper/2
      }
    }
    w.vec=as.vector(t(w.mat))
    as.vector(A%*%w.vec)
  }

  p.CG=p_func(t.vec, d.vec, group)

  ##### Jackknife #####
  N=length(t.vec)
  del=matrix(0,N,d)
  for(i in 1:N){ del[i,]=p_func(t.vec[-i],d.vec[-i],group[-i]) }
  jvar=(N-1)^2/N*cov(del)

  if(sum(is.na(jvar))>0){ P="NaN";Lower="NaN";Upper="NaN" }
  else{
    P=1-pnorm((p.CG-0.5)/sqrt(diag(jvar)), 0, 1)
    SE=sqrt(diag(jvar))
    Lower=pmax(0,p.CG-qnorm(1-0.05/2)*SE)
    Upper=pmin(p.CG-qnorm(0.05/2)*SE,1)
  }

  if(S.plot==TRUE){
    plot(survfit(Surv(t.vec,d.vec)~group),
         col=1:d,mark.time=mark.time,lwd=2,
         xlab="Time",ylab="Survival probability")
    abline(v=t.upper,lty="dotted")
    legend("topright",legend=1:d,lty=1,lwd=rep(2,d),
           col=1:d,bty="n")
    text(0.95*t.upper,0.01,expression(tau))
  }

  ### F-test ###
  if(is.null(C)){  C=diag(rep(1,d))-rep(1,d)%*%t(rep(1,d))/d  }
  T.mat=t(C)%*%ginv(C%*%t(C))%*%C
  TV=N*T.mat%*%jvar
  F.test=as.vector(N*t(p.CG)%*%T.mat%*%p.CG/sum(diag(TV)))

  ### Critical values ###
  Chi.mat=matrix(rchisq(R*d,df=1),R,d)
  Q.simu=Chi.mat%*%(Re(eigen(TV)$values)/sum(diag(TV)))
  c.simu=c("0.10"=sort(Q.simu)[(1-0.10)*R],
         "0.05"=sort(Q.simu)[(1-0.05)*R],
         "0.01"=sort(Q.simu)[(1-0.01)*R])

  df=sum(diag(TV))^2/sum(diag(TV%*%TV))
  c.anal=c("0.10"=qchisq(1-0.10,df=df)/df,
        "0.05"=qchisq(1-0.05,df=df)/df,
        "0.01"=qchisq(1-0.01,df=df)/df)

  P.value=c("simu"=mean(Q.simu>F.test),
            "anal"=1-pchisq(F.test*df,df=df))

  p.res=cbind(estimate=p.CG,se=SE,lower=Lower,upper=Upper,P)
  list(copula.parameter=alpha,p=p.res,Var=jvar,F=F.test,
       c.simu=c.simu,c.anal=c.anal,P.value=P.value)
}
