## ----eval=TRUE----------------------------------------------------------------
#install.packages("invgamstochvol")

## ----eval=TRUE----------------------------------------------------------------
library(invgamstochvol)

## ----eval=TRUE----------------------------------------------------------------
##simulate data
n=150
dat<-data.frame(Ydep=runif(n,0.3,1.4))
Ydep <- as.matrix(dat, -1,ncol=ncol(dat))
littlerho=0.95
r0=1
rho=diag(r0)*littlerho
p=4
n=4.1
T=nrow(Ydep)
Xdep <- Ydep[p:(T-1),]
if (p>1){
  for(lagi in 2:p){
    Xdep <- cbind(Xdep, Ydep[(p-lagi+1):(T-lagi),])
  }
}
T=nrow(Ydep)
Ydep <- as.matrix(Ydep[(p+1):T,])
T=nrow(Ydep)
unos <- rep(1,T)
Xdep <- cbind(unos, Xdep)

## ----eval=TRUE----------------------------------------------------------------
## obtain residuals
bOLS <- solve(t(Xdep) %*% Xdep) %*% t(Xdep) %*% Ydep
Res= Ydep- Xdep %*% bOLS
Res=Res[1:T,1]
b2=solve(t(Res) %*% Res/T) %*% (1-rho %*% rho)/(n-2)
Res=as.matrix(Res,ncol=1)

## ----eval=TRUE----------------------------------------------------------------
## obtain the log likelihood
LL1=lik_clo(Res,b2,n,rho)
LL1[1]

## ----eval = TRUE--------------------------------------------------------------
##Example using US data 
data1 <- US_Inf_Data
Ydep <- as.matrix(data1)
littlerho=0.95
r0=1
rho=diag(r0)*littlerho
p=4
n=4.1
T=nrow(Ydep)
Xdep <- Ydep[p:(T-1),]
if (p>1){
  for(lagi in 2:p){
    Xdep <- cbind(Xdep, Ydep[(p-lagi+1):(T-lagi),])
  }
}
T=nrow(Ydep)
Ydep <- as.matrix(Ydep[(p+1):T,])
T=nrow(Ydep)
unos <- rep(1,T)
Xdep <- cbind(unos, Xdep)

## ----eval=TRUE----------------------------------------------------------------
## obtain residuals
bOLS <- solve(t(Xdep) %*% Xdep) %*% t(Xdep) %*% Ydep
Res= Ydep- Xdep %*% bOLS
Res=Res[1:T,1]
b2=solve(t(Res) %*% Res/T) %*% (1-rho %*% rho)/(n-2)
Res=as.matrix(Res,ncol=1)

## ----eval=TRUE----------------------------------------------------------------
##obtain the log likelihood 
LL1=lik_clo(Res,b2,n,rho)
LL1[1]

## ----eval=TRUE----------------------------------------------------------------
deg=200
niter=200
AllSt=matrix(unlist(LL1[3]), ncol=1)
allctil=matrix(unlist(LL1[4]),nrow=T, ncol=(deg+1))
donde=(niter>deg)*niter+(deg>=niter)*deg 
alogfac=matrix(unlist(LL1[5]),nrow=(deg+1),ncol=(donde+1))
alogfac2=matrix(unlist(LL1[6]), ncol=1)
alfac=matrix(unlist(LL1[7]), ncol=1)

## ----eval=TRUE----------------------------------------------------------------
milaK=0
repli=5
keep0=matrix(0,nrow=repli, ncol=1)
for (jj in 1:repli)
{
  laK=DrawK0(AllSt,allctil,alogfac, alogfac2, alfac, n, rho, b2,nproc2=2)
  
  milaK=milaK+1/laK*(1/repli)
  keep0[jj]=mean(1/laK)/b2
}
ccc=1/b2
fefo=as.vector(milaK)*ccc

## ----eval=TRUE----------------------------------------------------------------
##obtain moving average of squared residuals
mRes=matrix(0,nrow=T,ncol=1)
Res2=Res*Res
bandi=5
for (iter in 1:T)
{  low=(iter-bandi)*(iter>bandi)+1*(iter<=bandi)
   up=(iter+bandi)*(iter<=(T-bandi))+T*(iter>(T-bandi))
  mRes[iter]=mean(Res2[low:up])
}

##plot the results
plot(fefo,type="l", col = "red", xlab="Time",ylab="Volatility Means")
lines(mRes, type="l", col = "blue")
legend("topright", legend = c("Stochastic Volatility", "Squared Residuals"),
       col = c("red", "blue"), lty = 1, cex = 0.8)


##usage of ourgeo to evaluate a 2F1 hypergeometric function
ourgeo(1.5,1.9,1.2,0.7)

