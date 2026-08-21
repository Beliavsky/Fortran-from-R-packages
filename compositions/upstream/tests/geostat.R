# source("acheck.R")

# try(library(compositions,lib.loc="compositions.Rcheck"))
# options(error=dump.frames)
options(warn=1)
library(compositions)
data(juraset)

X <- with(juraset,cbind(X,Y))
comp <- acomp(juraset,c("Cd","Cu","Pb","Co","Cr"))


lrv <- logratioVariogram(comp,X,maxdist=1,nbins=10)

plot(lrv)




fff <- CompLinModCoReg(~nugget()+sph(0.5)+exp(0.7),comp)
#fff <- CompLinModCoReg(~nugget()+sph(0.5),comp)
fit <- fit.lmc(lrv,fff,iterlim=100,print.level=0)
fit
fff(1:3)
plot(lrv,lrvg=vgram2lrvgram(fit))

x <- (0:60/60)*6
y <- (0:40/40)*6
Xnew <- cbind(rep(x,length(y)),rep(y,each=length(x)))


# Grid example
erg <- compOKriging(comp,X,Xnew,fit,err=FALSE)
image(x,y,structure(unclass(balance(erg$Z,~Cd/Cu)),dim=c(length(x),length(y))), asp=1)
par(mar=c(0,0,1,0))
suppressWarnings(
  pairwisePlot(erg$Z,
               panel=function(a,b,xlab,ylab){
                 image(x,y, asp=1, structure(log(a/b),dim=c(length(x),length(y))),main=paste("log(",xlab,"/",ylab,")",sep="")
                       )
                 points(X,pch=".")
                 }
               )
)

#erg <- compOKriging(comp,X,Xnew,fit$vg,err=TRUE)


# Checking 
ergR <- compOKriging(comp,X,X+1E-7,fit,err=FALSE)
plot(ergR$Z-comp) 
ergR <- compOKriging(comp,X,X,fit,err=FALSE)
plot(ergR$Z-comp) 

pairwisePlot(ilr(comp),ilr(ergR$Z))

ergR <- compOKriging(comp,X,X[rev(1:31),],fit,err=FALSE)
pairwisePlot(ilr(comp)[rev(1:31),],ilr(ergR$Z))

plot(ergR$Z)
plot(comp)
#

                                        # Small example
X <- rbind(c(0,1),c(1,0),c(1,1),c(0,0))
Z <- acomp(rbind(c(2,1,1),c(1,2,1),c(1,1,2),c(1,1,1)))
Xnew <- rbind(c(0,1),c(0.5,0.5))

vg <- CompLinModCoReg(~nugget()+sph(3),Z)
vg
(ergR <- compOKriging(Z,X,X,vg,err=FALSE))

## Explicit vgmodel; check the following functions prefixed with compositions::: 
#vgmodel <- function(h,nugget=0.2*parameterPosDefMat(diag(5)),sill=0.6*parameterPosDefMat(diag(5)),range1=1) {
#  (h>0) %o% parametricPosdefMat(nugget)+ vgram.sph(h,range=range1)%o%parametricPosdefMat(sill)
#}
# plot(lrv,lrvg=vgram2lrvgram(vgmodel))

#fit <- vgmFit(lrv,vgmodel)
#fit <- vgmFit(lrv,vgmodel,mode="ls")
#fit
#vgmGof(emp=lrv,vg=vgmodel,mode="ls")
#vgmGof(emp=lrv,vg=vgmodel,mode="log")
#lrv$h
#lrv$vg
#vgmGof(emp=lrv,vg=fit$vg,mode="ls")
#vgmGof(emp=lrv,vg=fit$vg,mode="log")

#plot(lrv,lrvg=vgram2lrvgram(fit$vg))



