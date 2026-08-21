\name{as.data.frame}
\alias{as.data.frame.acomp}
\alias{as.data.frame.rcomp}
\alias{as.data.frame.aplus}
\alias{as.data.frame.rplus}
\alias{as.data.frame.rmult}
\alias{as.data.frame.ccomp}
\alias{as.matrix.rmult}
\title{Convert "compositions" classes to data frames or matrices}
\description{
  Convert a compositional object to a dataframe
}
\usage{
\method{as.data.frame}{acomp}(x,\dots)
\method{as.data.frame}{rcomp}(x,\dots)
\method{as.data.frame}{aplus}(x,\dots)
\method{as.data.frame}{rplus}(x,\dots)
\method{as.data.frame}{rmult}(x,\dots)
\method{as.data.frame}{ccomp}(x,\dots)
\method{as.matrix}{rmult}(x,\dots)

          }
\arguments{
  \item{x}{an object to be converted to a dataframe}
  \item{\dots}{additional arguments are not used}
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\value{
  a data frame containing the given data, or (for rmult only) as matrix.
}
\examples{
data(SimulatedAmounts)
as.data.frame(acomp(sa.groups))
# The central perpose of providing this command is that the following
# works properly:
data.frame(acomp(sa.groups),groups=sa.groups.area)
}
\keyword{multivariate}
