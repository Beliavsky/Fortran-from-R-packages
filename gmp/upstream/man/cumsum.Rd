\name{cumsum}
\alias{cumsum.bigz}
\alias{cumsum.bigq}
\alias{sum.bigz}
\alias{sum.bigq}
\alias{prod.bigz}
\alias{prod.bigq}

\title{(Cumulative) Sums, Products of Large Integers and Rationals}
\description{
  Theses are methods to \sQuote{overload} the \code{\link{sum}()},
  \code{\link{cumsum}()} and \code{\link{prod}()} functions for big
  rationals and big integers.
}
\usage{
\S3method{cumsum}{bigz}(x)
\S3method{cumsum}{bigq}(x)
\S3method{sum}{bigz}(..., na.rm = FALSE)
\S3method{sum}{bigq}(..., na.rm = FALSE)
\S3method{prod}{bigz}(..., na.rm = FALSE)
\S3method{prod}{bigq}(..., na.rm = FALSE)
}
\arguments{
  \item{x, ...}{\R objects of class \code{bigz} or \code{bigq} or
    \sQuote{simple} numbers.}
  \item{na.rm}{logical indicating if missing values (\code{\link{NA}})
    should be removed before the computation.}
}
\value{
  return an element of class bigz or bigq.
}
\author{Antoine Lucas}
\seealso{\code{\link{apply}}
}

\examples{
 x <- as.bigz(1:12)
 cumsum(x)
 prod(x)
 sum(x)

 x <- as.bigq(1:12)
 cumsum(x)
 prod(x)
 sum(x)
}

\keyword{arith}

