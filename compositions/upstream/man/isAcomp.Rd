\name{is.acomp}
\alias{is.acomp}
\alias{is.aplus}
\alias{is.ccomp}
\alias{is.rcomp}
\alias{is.rplus}
\alias{is.rmult}
\title{Check for compositional data type}
\description{
  \code{is.}\var{XXX} returns \code{TRUE} if and only if its argument is
  of type \var{XXX} }
\usage{
is.acomp(x)
is.rcomp(x)
is.aplus(x)
is.rplus(x)
is.rmult(x)
is.ccomp(x)
          }
\arguments{
  \item{x}{any object to be checked}
}
\value{
  TRUE or FALSE
}
\details{
  These functions only check for the class of the object.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{acomp}},  \code{\link{rcomp}}
  \code{\link{aplus}},  \code{\link{rplus}}

}
\examples{
is.acomp(1:3)
is.acomp(acomp(1:3))
is.rcomp(acomp(1:3))
is.acomp(acomp(1:3)+acomp(1:3))
}
\keyword{classes}
