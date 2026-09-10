\name{frexpZ}
\title{Split Number into Fractional and Exponent of 2 Parts}
\alias{frexpZ}
\alias{frexp}% for search -- FIXME? make frexp() a generic with bigz() method?
\description{
  Breaks the number \code{x} into its binary significand
  (\dQuote{fraction}) \eqn{d \in [0.5, 1)}{d in [0.5, 1)}
  and \eqn{ex}, the integral exponent for 2, such that
  \eqn{x = d \cdot 2^{ex}.}{x = d * 2 ^ ex.}

  If \code{x} is zero, both parts (significand and exponent) are zero.
}
\usage{
frexpZ(x)
}
\arguments{
  \item{x}{integer or big integer (\code{\link{bigz}}).}
}
\value{
  a \code{\link{list}} with the two components
  \item{d}{a numeric vector whose absolute values are either zero, or in
    \eqn{[\frac{1}{2}, 1)}{[0.5, 1)}.}
  \item{exp}{an integer vector of the same length;
    note that \code{exp == 1 + floor(log2(x))}, and hence always
    \code{exp > log2(x)}.}
}
\author{Martin Maechler}
\seealso{
  \code{\link{log2}}, etc; for \code{\link{bigz}} objects built on
  (the C++ equivalent of) \code{frexp()}, actually GMP's
  \samp{mpz_get_d_2exp()}.
}
\examples{
frexpZ(1:10)
## and confirm :
with(frexpZ(1:10),  d * 2^exp)
x <- rpois(1000, lambda=100) * (1 + rpois(1000, lambda=16))
X <- as.bigz(x)
stopifnot(all.equal(x, with(frexpZ(x), d* 2^exp)),
          1+floor(log2(x)) == (fx <- frexpZ(x)$exp),
          fx == frexpZ(X)$exp,
          1+floor(log2(X)) == fx
)
}
\keyword{arith}
