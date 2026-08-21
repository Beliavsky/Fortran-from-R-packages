\name{norm}
\alias{norm}
\alias{norm.default}
\alias{norm.acomp}
\alias{norm.aplus}
\alias{norm.rcomp}
\alias{norm.rplus}
\alias{norm.rmult}
\alias{norm.matrix}
\title{Vector space norm}
\description{
  Each of the considered space structures has an associated norm, which is
  computed for each element by these functions.
}
\usage{
\method{norm}{default}(x,\dots)
\method{norm}{acomp}(x,\dots)
\method{norm}{rcomp}(x,\dots)
\method{norm}{aplus}(x,\dots)
\method{norm}{rplus}(x,\dots)
\method{norm}{rmult}(x,\dots)
\method{norm}{rmult}(x,\dots)
          }
\arguments{
  \item{x}{a dataset or a single vector of some type}
  \item{\dots}{currently not used, intended to select a 
              different norm rule in the future}
}
\value{
  The norms of the given vectors.
  ATTENTON: \code{norm.matrix} is a wrapper around base::norm
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{normalize}}
}
\examples{
data(SimulatedAmounts)
tmp <- acomp(sa.lognormals)
mvar(tmp)
sum(norm( tmp - mean(tmp) )^2)/(nrow(tmp)-1)

}
\keyword{multivariate}
