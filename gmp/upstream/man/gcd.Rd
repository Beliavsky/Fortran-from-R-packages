\name{gcd.bigz}
\title{Greatest Common Divisor (GCD) and Least Common Multiple (LCM)}
\alias{gcd.bigz}
\alias{gcd}
\alias{gcd.default}
\alias{lcm.default}
%\alias{lcm}% FIXME (R !):  conflict with graphics::lcm  is the trivial function
\alias{lcm.bigz}
\concept{GCD}
\concept{LCM}
\description{
  Compute the greatest common divisor (GCD) and least common multiple
  (LCM) of two (big) integers.
}
\usage{
\method{gcd}{bigz}(a, b)
lcm.bigz(a, b)
}
\arguments{
  \item{a,b}{Either integer, numeric, \code{\link{bigz}} or a string
    value; if a string, either starting with \code{0x} for hexadecimal,
    \code{0b} for binary or without prefix for decimal values.}
}
\value{
  An element of class bigz
}
\references{The GNU MP Library, see \url{https://gmplib.org}}
\author{Antoine Lucas}

\seealso{\code{\link{gcdex}}}

\examples{
gcd.bigz(210,342) # or also
lcm.bigz(210,342)
a <- 210 ; b <- 342
stopifnot(gcd.bigz(a,b) * lcm.bigz(a,b) == a * b)

## or
(a <- as.bigz("82696155787249022588"))
(b <- as.bigz("65175989479756205392"))
gcd(a,b) # 4
stopifnot(gcd(a,b) * lcm.bigz(a,b) == a * b)
}
\keyword{arith}

