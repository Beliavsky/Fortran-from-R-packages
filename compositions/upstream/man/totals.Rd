\name{totals}
\alias{totals}
\alias{totals.acomp}
\alias{totals.aplus}
\alias{totals.rcomp}
\alias{totals.rplus}
\alias{totals.ccomp}
\title{Total sum of amounts }
\description{
  Calculates the total amount by summing the individual parts. 
}
\usage{
totals(x,\dots)
\method{totals}{acomp}(x,\dots,missing.ok=TRUE)
\method{totals}{rcomp}(x,\dots,missing.ok=TRUE)
\method{totals}{aplus}(x,\dots,missing.ok=TRUE)
\method{totals}{rplus}(x,\dots,missing.ok=TRUE)
\method{totals}{ccomp}(x,\dots,missing.ok=TRUE)
          }
\arguments{
  \item{x}{an amount/amount dataset}
  \item{\dots}{not used, only here for generic purposes}
  \item{missing.ok}{if TRUE ignores missings; if FALSE issues an error if
    the total cannot be calculated due to missings.}
}
\section{Missing Policy}{
    if \code{missing.ok=TRUE} missings are just regarded as 0, if
    \code{missing.ok=FALSE} WZERO values is still regarded as 0 and other
    sorts lead to \code{NA} in the respective totals.
}
\value{
  a numeric vector of length equal to \code{ncol(x)} containing the total amounts
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{aplus}}
}
\examples{
data(SimulatedAmounts)
totals(acomp(sa.lognormals))
totals(rcomp(sa.lognormals,total=100))
totals(aplus(sa.lognormals))
totals(rplus(sa.lognormals))
aplus(acomp(sa.lognormals),total=totals(aplus(sa.lognormals)))


}
\keyword{multivariate}
