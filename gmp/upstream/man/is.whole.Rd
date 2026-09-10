\name{is.whole}
\alias{is.whole}
\alias{is.whole.default}
\alias{is.whole.bigz}
\alias{is.whole.bigq}
\title{Whole ("Integer") Numbers}
\description{
  Check which elements of \code{x[]} are integer valued aka
  \dQuote{whole} numbers.
}
\usage{
is.whole(x)
\S3method{is.whole}{default}(x)
\S3method{is.whole}{bigz}(x)
\S3method{is.whole}{bigq}(x)
}
\arguments{
  \item{x}{any \R vector}
}
\value{
  logical vector of the same length as \code{x}, indicating where
  \code{x[.]} is integer valued.
}
\author{Martin Maechler}
\seealso{
 \code{\link{is.integer}(x)} (\pkg{base} package) checks for the
 \emph{internal} mode or class; not if \code{x[i]} are integer valued.

 The \code{\link[Rmpfr]{is.whole}()} method for "mpfr" numbers.
}
\examples{
 is.integer(3) # FALSE, it's internally a double
 is.whole(3)   # TRUE
 ## integer valued complex numbers  (two FALSE) :
 is.whole(c(7, 1 + 1i, 1.2, 3.4i, 7i))
 is.whole(factorialZ(20)^(10:12)) ## "bigz" are *always* whole numbers
 q <- c(as.bigz(36)^50 / as.bigz(30)^40, 3, factorialZ(30:31), 12.25)
 is.whole(q) # F T T T F
}
\keyword{math}
