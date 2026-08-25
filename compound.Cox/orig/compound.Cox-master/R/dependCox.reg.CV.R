dependCox.reg.CV <-
function(t.vec,d.vec,X.mat,K=5,G=20){

X.mat=as.matrix(X.mat)
### Grid search on CV ####
tau_grid=seq(0.0001,0.9,length=G)
alpha_grid=2*tau_grid/(1-tau_grid)
C_grid=NULL
for(i in 1:G){
  C_grid[i]=cindex.CV(t.vec,d.vec,X.mat,alpha_grid[i])
}
alpha=alpha_grid[C_grid==max(C_grid)][1]

######### univariate Cox with dependent censoring #########
p=ncol(X.mat)
Beta=SE=Z=P=numeric(p)

for(j in 1:p){
  res=dependCox.reg(t.vec,d.vec,X.mat[,j],alpha=alpha,var=TRUE)
  Beta[j]=res[1]
  SE[j]=res[2]
}
Z=Beta/SE
P=1-pchisq(Z^2,df=1)

plot(tau_grid,C_grid,xlab="Kendall's tau ( alpha )",ylab="CV( alpha )",type="b",lwd=3)
points(tau_grid[C_grid==max(C_grid)][1],max(C_grid),col="red",pch=17,cex=2)

names(Z)=names(P)=names(Beta)=names(SE)=colnames(X.mat)
list(beta=Beta,SE=SE,Z=Z,P=P,alpha=alpha,c_index=max(C_grid))

}
