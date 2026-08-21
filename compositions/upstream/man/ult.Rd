\name{ult}
\alias{ult}
\alias{ultInv}
\alias{Kappa}
\title{Uncentered log transform}
\description{
  Compute the uncentered log ratio transform of a (dataset of)
  composition(s) and its inverse.
}
\usage{
          ult( x ,...)
          ultInv( z ,..., orig=gsi.orig(z))
          Kappa( x ,...)
          }
\arguments{
  \item{x}{a composition or a data matrix of compositions, not necessarily closed}
  \item{z}{the ult-transform of a composition or
    clr-transforms of compositions (or a data matrix), not necessarily centered}
  \item{\dots}{for generic use only}
  \item{orig}{a compositional object which should be mimicked
             by the inverse transformation. It is the generic
	     argument. Typically the \code{orig} argument is stored as an attribute 
	     in \code{z} and will be extracted automatically by this function; if this
	     fails, \code{orig} can be set equal to the dataset that
	     was transformed in the first place.}
}
\value{
  \code{ult} gives the uncentered log transform,\cr
  \code{ultInv} gives closed compositions with the given
  ult/clr-transforms\cr
  \code{Kappa} gives the difference between the clr and the ult
  transforms. It is quite linked to information measures.
  }
\details{
  The ult-transform is simply the elementwise log of the closed
  composition. The ult has some important properties in the scope 
  of Information Theory of probability vectors (but might be mostly 
  misleading for exploratory analysis of compositions). DO NOT USE if 
  you do not know what you are doing.

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{ilr}},\code{\link{alr}},\code{\link{apt}}
}
\examples{
(tmp <- ult(c(1,2,3)))
ultInv(tmp)
ultInv(tmp) - clo(c(1,2,3)) # 0
data(Hydrochem)
cdata <- Hydrochem[,6:19]
pairs(ult(cdata),pch=".")
Kappa(c(1,2,3))
}
\keyword{multivariate}
