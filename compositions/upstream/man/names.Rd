\name{names}
\alias{names.acomp}
\alias{names.aplus}
\alias{names.rcomp}
\alias{names.rplus}
\alias{names.rmult}
\alias{names.ccomp}
\alias{names<-.acomp}
\alias{names<-.ccomp}
\alias{names<-.aplus}
\alias{names<-.rcomp}
\alias{names<-.rplus}
\alias{names<-.rmult}
\title{The names of the parts }
\description{
  The \code{names} function provide a transparent way to access the names of
  the parts regardless of the shape of the dataset or data item.
}
\usage{
\method{names}{acomp}(x)
\method{names}{rcomp}(x)
\method{names}{aplus}(x)
\method{names}{rplus}(x)
\method{names}{rmult}(x)
\method{names}{ccomp}(x)
\method{names}{acomp}(x) <- value
\method{names}{rcomp}(x) <- value
\method{names}{aplus}(x) <- value
\method{names}{rplus}(x) <- value
\method{names}{rmult}(x) <- value
\method{names}{ccomp}(x) <- value

          }
\arguments{
  \item{x}{an amount/amount dataset}
  \item{value}{the new names of the parts}
  %\item{\dots}{not used, only here for generics}
}
\value{
  a character vector giving the names of the parts
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{aplus}}
}
\examples{
data(SimulatedAmounts)
tmp <- acomp(sa.lognormals)
names(tmp)
names(tmp) <- c("x","y","z")
tmp
}
\keyword{multivariate}
