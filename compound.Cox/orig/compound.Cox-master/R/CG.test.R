CG.test <-
  function(t.vec,d.vec,PI,cutoff=median(PI),
           alpha=2,copula=CG.Clayton,
           S.plot=TRUE,N=10000,mark.time=TRUE){

### Estimating survival for two groups ###
CG.surv=function(t.vec,d.vec,PI,
                 cutoff=median(PI),Plot=S.plot){
  group=(PI<=cutoff) ## Good or poor group

  t.good=t.vec[group]
  d.good=d.vec[group]
  res.good=copula(t.good,d.good,alpha,S.plot=FALSE)
  S.good=res.good$surv

  t.poor=t.vec[!group]
  d.poor=d.vec[!group]
  res.poor=copula(t.poor,d.poor,alpha,S.plot=FALSE)
  S.poor=res.poor$surv

  Tau=min(max(t.poor),max(t.good)) ## Pepe & Fleming (1989, 1991)
  t.Tau=sort(t.vec)[sort(t.vec)<=Tau]
  M=sum(t.vec<=Tau)
  group.o=group[order(t.vec)]

  SS.good=SS.poor=1
  for(i in 1:M){
    if(group.o[i]==TRUE){
      SS.good[i]=S.good[sum(group.o[1:i])]
      SS.poor[i]=SS.poor[max(1,i-1)]
    }else{
      SS.poor[i]=S.poor[sum((!group.o)[1:i])]
      SS.good[i]=SS.good[max(1,i-1)]
    }
  }

  RMS.good=sum( c(1,SS.good[1:(M-1)])*diff(c(0,t.Tau)) )
  RMS.poor=sum( c(1,SS.poor[1:(M-1)])*diff(c(0,t.Tau)) )
  L1=sum( c(0,abs(SS.good[1:(M-1)]-SS.poor[1:(M-1)]))*diff(c(0,t.Tau)) )

  ### Plot survival functions ###
  if(Plot==TRUE){
    plot(c(0,sort(t.good)),c(1,S.good),type="s",lwd=2,
         xlim=c(0,max(t.vec)),ylim=c(0,1),
         col="blue",ylab="Survival probability",xlab="Time")
    points(c(0,sort(t.poor)),c(1,S.poor),type="s",lwd=2,col="red")
    if(mark.time==TRUE){
      points(sort(t.good[d.good==0]),S.good[d.good[order(t.good)]==0],
             pch=3,cex=1,col="blue")
      points(sort(t.poor[d.poor==0]),S.poor[d.poor[order(t.poor)]==0],
             pch=3,cex=1,col="red")
    }

    polygon(c(t.Tau,rev(t.Tau)),c(SS.poor,rev(SS.good)),
            density=40,col="thistle1",border=NA)

    legend("topright",legend="Good prognosis (PI <= c)",col="blue",
           lty=1,lwd=2,pch=3,bg="transparent",bty="n")
    legend("bottomleft",legend="Poor prognosis (PI > c)",col="red",
           lty=1,lwd=2,pch=3,bg="transparent",bty="n")
    legend("bottomright",legend=c("Kendall's tau = ",round(res.good$tau,2)),
           bg="transparent",bty="n")
    abline(v=Tau,col=rgb(0,0,0, alpha=0.1),cex=2)

  }

  list(Tau=Tau,RMS.poor=RMS.poor,L1=L1,t.poor=res.poor$time,
       d.poor=d.poor[order(t.poor)],S.poor=S.poor,
       RMS.good=RMS.good,t.good=res.good$time,
       d.good=d.good[order(t.good)],S.good=S.good)
}

fit=CG.surv(t.vec,d.vec,PI,cutoff=cutoff)
RMSD=fit$RMS.good-fit$RMS.poor
MD=RMSD/fit$Tau
IntegratedL1=fit$L1
ML1=fit$L1/fit$Tau

#### Permutation test #####
n=length(PI)
MD.perm=ML1.perm=numeric(N)
set.seed(1)

for(j in 1:N){
  perm=sample(size=n,1:n,replace=FALSE) #permutation
  fit.perm=CG.surv(t.vec,d.vec,PI[perm],
                   cutoff=cutoff,Plot=FALSE)
  RMSD.perm=fit.perm$RMS.good-fit.perm$RMS.poor
  MD.perm[j]=RMSD.perm/fit.perm$Tau
  ML1.perm[j]=fit.perm$L1/fit.perm$Tau
}
P_MD=mean( abs(MD.perm)>abs(MD) )
P_L1=mean( ML1.perm>ML1 )

if(S.plot==TRUE){
  legend("center",legend=c("P-value = ",P_MD),
       bg="transparent",bty="n")
}

list(test=c(Survival.diff=MD,RMSTD=fit$Tau*MD,P.value=P_MD),
     L1.test=c(L1.distance=ML1,Integrated.L1=IntegratedL1,P.value=P_L1),
     Good=c(n=sum(PI<=cutoff),n.event=sum(d.vec[PI<=cutoff]),
            RMST=fit$RMS.good,meanPI=mean(PI[PI<=cutoff])),
     Poor=c(n=sum(PI>cutoff),n.event=sum(d.vec[PI>cutoff]),
            RMST=fit$RMS.poor,meanPI=mean(PI[PI>cutoff])))
}


