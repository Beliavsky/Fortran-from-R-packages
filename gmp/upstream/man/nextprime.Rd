\name{nextprime}
\alias{nextprime}
\title{Next Prime Number}
\description{
  Return the next prime number, say \eqn{p}, with \eqn{p > n}{p > n}.
}
\usage{
nextprime(n)
}
\arguments{
  \item{n}{Integer}
}
\details{
  This function uses probabilistic algorithm to identify
  primes.  For practical purposes, it is adequate, the chance
  of a composite passing will be extremely small.
}
\value{
  A (probably) prime number
}
\references{The GNU MP Library, see \url{https://gmplib.org}}

\author{Antoine Lucas}

\seealso{
  \code{\link{isprime}} and its references and examples.
}
\examples{
nextprime(14)
## still very fast:
(p <- nextprime(1e7))
## to be really sure { isprime()  gives "probably prime" } :
stopifnot(identical(p, factorize(p)))
}
\keyword{arith}


