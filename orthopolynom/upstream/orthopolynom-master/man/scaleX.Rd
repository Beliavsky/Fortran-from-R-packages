\name{scaleX}
\alias{scaleX}
\title{ Scale values from [a,b] to [u.v] }
\description{
  This function returns a vector of values that have been mapped from the interval [a,b]
  to the interval [u.v].
}
\usage{
scaleX(x, a = min(x, na.rm = TRUE), b = max(x, na.rm = TRUE), u, v)
}
\arguments{
  \item{x}{ A numerical vector of values to be mapped into a target interval }
  \item{a}{ A numerical lower bound for the domain interval with min(x) as the default value }
  \item{b}{ A numerical upper bound for the domain interval with max(x) as the default value }
  \item{u}{ A numerical lower bound for the target interval }
  \item{v}{ A numerical upper bound for the target interval }
}
\details{
  Target lower and/or upper bounds can be \eqn{-\infty} and \eqn{\infty}, respectively. This accomodates
  finite target intervals, semi-infinite target intervals and infinite target intervals.
}
\value{
  A vector of transformed values with four attributes.  The first attribute is called "a" and
  it is the domain interval lower bound.  The second attribute is called "b" and
  it is the domain interval upper bound. The third attribute is called "u" and
  it is the target interval lower bound.  The fourth attribute is called
  "v" and it is the target interval upper bound.
}
\references{
Seber, G. A. F. (1997) Linear Regression Analysis, New York.
}
\author{
Frederick Novomestky \email{fnovomes@poly.edu},
Gregor Gorjanc \email{gregor.gorjanc@bfro-uni-lj.si}
}
\examples{
x <- rnorm( 1000, 0, 10 )
y0 <- scaleX( x, u=0 , v=1 )
y1 <- scaleX( x, u=-1, v=1 )
y2 <- scaleX( x, u=-Inf, v=0 )
y3 <- scaleX( x, u=0, v=Inf )
y4 <- scaleX( x, u=-Inf, v=Inf )
}
\keyword{ math }
