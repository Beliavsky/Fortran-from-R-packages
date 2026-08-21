\name{compOKriging}
\alias{compOKriging}
\title{Compositional Ordinary Kriging}
\description{
  Geostatistical prediction for compositional data with missing values.
}
\usage{
compOKriging(comp,X,Xnew,vg,err=FALSE)
          }
\arguments{
  \item{comp}{an acomp compositional dataset}
  \item{X}{A dataset of locations}
  \item{Xnew}{The locations, where a geostatistical prediction should be
    computed.}
  \item{vg}{A compositional variogram function.}
  \item{err}{boolean: If true kriging errors are computed additionally. A bug was found here; the argument is currently disabled.}
}
\value{
  A list of class \code{"logratioVariogram"}
  \item{X}{The new locations as given by Xnew}
  \item{Z}{The predicted values as acomp compositions.}
  \item{err}{A bug has been found here. This output is currently disabled (An ncol(Z)xDxD array with the clr kriging errors.)}
}
\details{
   
  The function performes multivariate ordinary kriging of compositions
  based on transformes addapted to the missings in every case. The
  variogram is assumed to be a clr variogram.
  
}
\references{
  Pawlowsky-Glahn, Vera and Olea, Ricardo A. (2004) Geostatistical
  Analysis of Compositional Data, Oxford University Press, Studies in
  Mathematical Geology
  
  Tolosana (2008) ...
  
  Tolosana, van den Boogaart, Pawlowsky-Glahn (2009) Estimating and
  modeling variograms of compositional data with occasional missing
  variables in R, StatGis09
  
  
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{vgram2lrvgram}},
  \code{\link{CompLinModCoReg}},
  \code{\link{vgmFit}}
}
\examples{
\dontrun{
# Load data
data(juraset)
X <- with(juraset,cbind(X,Y))
comp <- acomp(juraset,c("Cd","Cu","Pb","Co","Cr"))
lrv <- logratioVariogram(comp,X,maxdist=1,nbins=10)
plot(lrv)

# Fit a variogram model
vgModel <- CompLinModCoReg(~nugget()+sph(0.5)+R1*exp(0.7),comp)
fit <- vgmFit2lrv(lrv,vgModel)
fit
plot(lrv,lrvg=vgram2lrvgram(fit$vg))

# Define A grid
x <- (0:10/10)*6
y <- (0:10/10)*6
Xnew <- cbind(rep(x,length(y)),rep(y,each=length(x)))

# Kriging
erg <- compOKriging(comp,X,Xnew,fit$vg,err=FALSE)
par(mar=c(0,0,1,0))
pairwisePlot(erg$Z,panel=function(a,b,xlab,ylab) {image(x,y,
                     structure(log(a/b),dim=c(length(x),length(y))),
                     main=paste("log(",xlab,"/",ylab,")",sep=""));points(X,pch=".")})


# Check interpolation properties 
ergR <- compOKriging(comp,X,X,fit$vg,err=FALSE)
pairwisePlot(ilr(comp),ilr(ergR$Z))
ergR <- compOKriging(comp,X,X+1E-7,fit$vg,err=FALSE)
pairwisePlot(ilr(comp),ilr(ergR$Z))
ergR <- compOKriging(comp,X,X[rev(1:31),],fit$vg,err=FALSE)
pairwisePlot(ilr(comp)[rev(1:31),],ilr(ergR$Z))
}
}
\keyword{multivariate}
