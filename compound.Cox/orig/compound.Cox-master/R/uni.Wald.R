uni.Wald=function(t.vec, d.vec, X.mat){
  n=length(t.vec)
  X.mat=as.matrix(X.mat)
  p=ncol(X.mat)

  Beta=Z=P=numeric(p)
  for(j in 1:p){
    res=summary(coxph(Surv(t.vec,d.vec)~X.mat[,j]))$coef
    Beta[j]=res[1]
    Z[j]=res[4]
    P[j]=res[5]
  }
  names(Z)=names(P)=names(Beta)=colnames(X.mat)
  list(beta=Beta,Z=Z,P=P)
}
