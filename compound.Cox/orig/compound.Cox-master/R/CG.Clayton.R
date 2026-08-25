CG.Clayton <-
function(t.vec,d.vec,alpha,S.plot=TRUE,S.col="black"){

  alpha=max(alpha,0.000000001) ### prevent unstability for small alpha ###

  n=length(t.vec)
  R=n:1
  t.sort=sort(t.vec)
  d.sort=d.vec[order(t.vec)]
  A=( (R-1)/n )^(-alpha)-(R/n)^(-alpha)
  A[n]=0
  S=( 1+cumsum(A*d.sort) )^(-1/alpha)
  if(sum(S<=0.5)>=1){med=min(t.sort[S<=0.5])}else{med=Inf}

  if(S.plot==TRUE){
    plot(c(0,t.sort),c(1,S),type="s",lwd=3,
         ylim=c(0,1),ylab="Survival probability",xlab="Time",col=S.col)
    points(t.sort[d.sort==0],S[d.sort==0],pch=3,cex=2,col=S.col)
  }

  list( tau=alpha/(alpha+2), median=med, time=t.sort, n.risk=R, surv=S )
}
