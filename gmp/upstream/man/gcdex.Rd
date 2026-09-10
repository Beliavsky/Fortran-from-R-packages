\name{gcdex}
\alias{gcdex}

\title{Compute Bezoult Coefficient}
\description{
  Compute g,s,t as
  \eqn{as + bt = g = gcd(a,b)}{as + bt = g = gcd(a,b)}.
  s and t are also known as Bezoult coefficients.
}
\usage{
gcdex(a, b)
}

\arguments{
  \item{a,b}{either integer, numeric, character string, or of class \code{"bigz"};
    If a string, either starting with \code{"0x"} for hexadecimal, \code{"0b"} for
    binary or without prefix for decimal values.}
}

\value{
  a class \code{"bigz"} vector of length 3 with (long integer) values
  \eqn{g, s, t}.
}

\references{The GNU MP Library, see \url{https://gmplib.org}}
\author{Antoine Lucas}

\seealso{\code{\link{gcd.bigz}}}
\examples{
gcdex(342,654)
}
\keyword{arith}

