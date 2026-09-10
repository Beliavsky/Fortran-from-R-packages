\name{Bigq_operators}
\alias{add.bigq}
\alias{+.bigq}
\alias{sub.bigq}
\alias{-.bigq}
\alias{mul.bigq}
\alias{*.bigq}
\alias{div.bigq}
\alias{/.bigq}
\alias{pow.bigq}
\alias{^.bigq}
\alias{abs.bigq}

\title{Basic arithmetic operators for large rationals}
\description{
  Addition, subtraction, multiplication, division, and absolute value for
  large rationals, i.e. \code{"\link{bigq}"} class \R objects.
}
\usage{
add.bigq(e1, e2)
\method{+}{bigq}(e1, e2)

sub.bigq(e1, e2=NULL)
\method{-}{bigq}(e1, e2)

mul.bigq(e1, e2)
\method{*}{bigq}(e1, e2)

div.bigq(e1, e2)
\method{/}{bigq}(e1, e2)

\method{^}{bigq}(e1, e2)

\method{abs}{bigq}(x)
}
\arguments{
  \item{e1,e2, x}{of class \code{"\link{bigq}"}, or (\code{e1} and
    \code{e2}) integer or string from an integer}
}
\details{
  Operators can be use directly when the objects are of class \code{"\link{bigq}"}:
  a + b, a * b, etc, and \code{a ^ n}, where \code{n} must be coercable
  to a biginteger (\code{"\link{bigz}"}).
}
\value{
  A bigq class representing the result of the arithmetic operation.
}

\author{Immanuel Scholz and Antoine Lucas}

\examples{
## 1/3 + 1 = 4/3 :
as.bigq(1,3) + 1

r <- as.bigq(12, 47)
stopifnot(r ^ 3 == r*r*r)
}
\keyword{arith}

