\name{variograms}
\alias{vgram.sph}
\alias{vgram.exp}
\alias{vgram.gauss}
\alias{vgram.cardsin}
\alias{vgram.lin}
\alias{vgram.pow}
\alias{vgram.nugget}
\title{Variogram functions}
\description{
  Valid scalar variogram model functions.  
}
\usage{
vgram.sph( h , nugget = 0, sill = 1, range= 1,... )
vgram.exp( h , nugget = 0, sill = 1, range= 1,... )
vgram.gauss( h , nugget = 0, sill = 1, range= 1,... )
vgram.cardsin( h , nugget = 0, sill = 1, range= 1,... )
vgram.lin( h , nugget = 0, sill = 1, range= 1,... )
vgram.pow( h , nugget = 0, sill = 1, range= 1,... )
vgram.nugget( h , nugget = 1,...,tol=1E-8 )
}
\arguments{
  \item{h}{a vector providing distances, a matrix of distance vectors in
    its rows or a data.frame of distance vectors. }
  \item{nugget}{The size of the nugget effect (i.e. the limit to 0). At
    zero itself the value is always 0.}
  \item{sill}{The sill (i.e. the limit to infinity)}
  \item{range}{The range parameter. I.e. the distance in which sill is
    reached or if this does not exist, where the value is in some sense
    nearly the sill.}
  \item{\dots}{not used}
  \item{tol}{The distance that is considered as nonzero.}
}
\value{
  A vector of size NROW(h), giving the variogram values.
}
\details{
  The univariate variograms are used in the CompLinCoReg as building
  blocks of multivariate variogram models. 
  \describe{
    \item{sph}{Spherical variogram}
    \item{exp}{Exponential variogram}
    \item{gauss}{The Gaussian variogram.}
    \item{gauss}{The cardinal sine variogram.}
    \item{lin}{Linear Variogram. Increases over the sill, which is
      reached at \code{range}.}
    \item{pow}{The power variogram. Increases over the sill, which is
      reached at \code{range}. }
    \item{nugget}{The pure nugget effect variogram. }
  }
}
\references{
  Cressie, N.C. (1993) Spatial statistics
  
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
data(juraset)
X <- with(juraset,cbind(X,Y))
comp <- acomp(juraset,c("Cd","Cu","Pb","Co","Cr"))
lrv <- logratioVariogram(comp,X,maxdist=1,nbins=10)
plot(lrv)
}
}
\keyword{multivariate}
