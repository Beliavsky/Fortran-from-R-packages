CG.Frank <-
function(t.vec,d.vec,alpha,S.plot=TRUE,S.col="black"){

  if(abs(alpha)<0.000001){alpha=0.000001} ### prevent unstability for small alpha ###

  n=length(t.vec)
  R=n:1
  t.sort=sort(t.vec)
  d.sort=d.vec[order(t.vec)]
  A=log(  (exp(-alpha*(R-1)/n)-1)/(exp(-alpha*R/n)-1)  )
  A[n]=0
  S=-1/alpha*log(  1+(exp(-alpha)-1)*exp(cumsum(A* d.sort))  )
  if(sum(S<=0.5)>=1){med=min(t.sort[S<=0.5])}else{med=Inf}

  if(S.plot==TRUE){
    plot(c(0,t.sort),c(1,S),type="s",lwd=3,
         ylim=c(0,1),ylab="Survival probability",xlab="Time",col=S.col)
    points(t.sort[d.sort==0],S[d.sort==0],pch=3,cex=2,col=S.col)
  }
  func1=function(x){x/(exp(x)-1)}
  tau=1-4/alpha*(1-integrate(func1,0,alpha)$value/alpha)
  list( tau=tau, median=med, time=t.sort, n.risk=R, surv=S )
}
