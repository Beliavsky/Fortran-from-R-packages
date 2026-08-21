\name{lrvgram}
\alias{vgram2lrvgram}
\alias{cgram2vgram}
\title{vgram2lrvgram}
\description{
  Transforms model functions for different types of compositional (logratio)(co)variograms.
}
\usage{
cgram2vgram(cgram)
vgram2lrvgram(vgram)
}
\arguments{
  \item{cgram}{A (matrix valued) covariance function.}
  \item{vgram}{A (matrix valued) variogram functions.}
}
\value{
  A function that takes the same parameters as the input function
  (through a \dots parameterlist), but provides the correponding
  variogram values (cgram2vgram) or logratio Variogram (vgram2lrvgram) values. 
}
\details{
  The variogram is given by \code{cgram(0)-cgram(h)} and
  \code{lrvgram(h)[,i,j]==vgram(h)[,i,i]+vgram(h)[,i,j]-2*vgram(h)[,i,j]}.  

  The logratio-variogram is the set of variograms of each of the pairwise
  logratios. It can be proven that it carries the same information as a
  usual multivariate variogram. The great advantage is that all the
  funcitions have a direct interpreation and can be estimated even with
  (MAR) missings in the dataset. 
}
\references{
  Tolosana, van den Boogaart, Pawlowsky-Glahn (2009) Estimating and
  modeling variograms of compositional data with occasional missing
  variables in R, StatGis09
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{logratioVariogram}},
  \code{\link{CompLinModCoReg}},
  \code{\link{vgmFit}}
}
\examples{
data(juraset)
comp <- acomp(juraset,c("Cd","Cu","Pb","Co","Cr"))
vg <- CompLinModCoReg(~nugget()+sph(0.5)+R1*exp(0.7),comp)
vg(1:3)
vgram2lrvgram(vg)(1:3)
}
\keyword{multivariate}
