\name{normalize}
\alias{normalize}
\alias{normalize.default}
\title{Normalize vectors to norm 1}
\description{
  Normalize vectors to norm 1.
}
\usage{
normalize(x,\dots)
\method{normalize}{default}(x,\dots)
}
\arguments{
  \item{x}{a dataset or a single vector of some type}
  \item{\dots}{currently not used, intended to select 
               a different norm in the future}
}
\value{
  The vectors given, but normalized to norm 1.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{norm.rmult}}
}
\examples{
data(SimulatedAmounts)
normalize(c(1,2,3))
normalize(acomp(c(1,2,3)))
norm(normalize(acomp(sa.groups)))
}
\keyword{multivariate}
