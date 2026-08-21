\name{split}
\alias{split.acomp}
\alias{split.aplus}
\alias{split.rcomp}
\alias{split.rplus}
\alias{split.rmult}
\alias{split.ccomp}
\title{Splitting datasets in groups given by factors}
\description{
  Splits data sets of compositions in groups given by factors, 
  and gives the same class as the data to the result.
}
\usage{
\method{split}{acomp}(x,f,drop=FALSE,\dots)
\method{split}{rcomp}(x,f,drop=FALSE,\dots)
\method{split}{aplus}(x,f,drop=FALSE,\dots)
\method{split}{rplus}(x,f,drop=FALSE,\dots)
\method{split}{rmult}(x,f,drop=FALSE,\dots)
\method{split}{ccomp}(x,f,drop=FALSE,\dots)
          }
\arguments{
  \item{x}{a dataset or a single vector of some type}
  \item{f}{a factor that defines the grouping or a list of factors }
  \item{drop}{drop=FALSE also gives (empty) datsets for empty
    categories}
  \item{\dots}{Further arguments passed to split.default.
    Currently (and probably) without any use.}
}
\value{
  a list of objects of the same type as \code{x}.
}
\examples{
  data(SimulatedAmounts)
  split(acomp(sa.groups),sa.groups.area)
  lapply( split(acomp(sa.groups),sa.groups.area), mean)
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link[base]{split}}
}
\keyword{multivariate}
