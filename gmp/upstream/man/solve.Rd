\name{solve.bigz}
\alias{solve.bigz}
\alias{solve.bigq}

\title{Solve a system of equation}
\description{
  This generic function solves the equation \eqn{a \%*\% x = b} for
  \eqn{x}, where \eqn{b} can be either a vector or a matrix.

  If a and b are rational, return is a rational matrix.

  If a and b are big integers (of class bigz) solution is in Z/nZ if
  there is a common modulus, or a rational matrix if not.
}
\usage{
\method{solve}{bigz}(a, b, ...)
\method{solve}{bigq}(a, b, ...)
}
\arguments{
  \item{a,b}{A element of class bigz or bigq}
  \item{\dots}{Unused}
}
\details{
  It uses the Gauss and trucmuch algo \dots (to be detailled).%% FIXME
}
\value{
  If a and b are rational, return is a rational matrix.

  If a and b are big integers (of class bigz) solution is in Z/nZ if
  there is a common modulus, of a rational matrix if not.
}

\author{Antoine Lucas}

\seealso{\code{\link{solve}}}

\examples{
x <- matrix(1:4,2,2)  ## standard solve :
solve(x)

q <- as.bigq(x) ## solve with rational
solve(q)

z <- as.bigz(x)
modulus(z) <- 7  ## solve in Z/7Z :
solve(z)

b <- c(1,3)
solve(q,b)
solve(z,b)

## Inversion of ("non-trivial") rational matrices :

A <- rbind(c(10, 1,  3),
           c( 4, 2, 10),
           c( 1, 8,  2))
(IA.q <- solve(as.bigq(A))) # fractions..
stopifnot(diag(3) == A \%*\% IA.q)# perfect

set.seed(5); B <- matrix(round(9*runif(5^2, -1,1)), 5)
B
(IB.q <- solve(as.bigq(B)))
stopifnot(diag(5) == B \%*\% IB.q, diag(5) == IB.q \%*\% B,
          identical(B, asNumeric(solve(IB.q))))
}
\keyword{arith}

