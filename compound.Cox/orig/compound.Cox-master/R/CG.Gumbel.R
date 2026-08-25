CG.Gumbel <-
function(t.vec,d.vec,alpha,S.plot=TRUE,S.col="black"){

  alpha=max(alpha,0) ### negative alpha is not allowed ###
  n=length(t.vec)
  R=n:1
  t.sort=sort(t.vec)
  d.sort=d.vec[order(t.vec)]
  A=( -log((R-1)/n) )^(alpha+1)-( -log(R/n) )^(alpha+1)
  A[n]=0
  S=exp( -( cumsum(A*d.sort) )^(1/(1+alpha)) )
  if(sum(S<=0.5)>=1){med=min(t.sort[S<=0.5])}else{med=Inf}

  if(S.plot==TRUE){
    plot(c(0,t.sort),c(1,S),type="s",lwd=3,
         ylim=c(0,1),ylab="Survival probability",xlab="Time",col=S.col)
    points(t.sort[d.sort==0],S[d.sort==0],pch=3,cex=2,col=S.col)
  }

  list( tau=(alpha/(alpha+1)), median=med, time=t.sort, n.risk=R, surv=S )
}
