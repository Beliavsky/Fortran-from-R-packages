\name{Extremes}
\title{Extrema (Maxima and Minima)}

\alias{max.bigz}
\alias{max.bigq}
\alias{min.bigz}
\alias{min.bigq}
\alias{which.max,bigq-method}
\alias{which.max,bigz-method}
\alias{which.min,bigq-method}
\alias{which.min,bigz-method}

\usage{
\method{max}{bigz}(..., na.rm=FALSE)
\method{max}{bigq}(..., na.rm=FALSE)
\method{min}{bigz}(..., na.rm=FALSE)
\method{min}{bigq}(..., na.rm=FALSE)

\S4method{which.min}{bigz}(x)
%S4method{which.min}{bigq}(x)
%S4method{which.max}{bigz}(x)
\S4method{which.max}{bigq}(x)
}

\description{
  We provide S3 \code{\link{methods}} for \code{\link{min}} and
  \code{\link{max}} for big rationals (\code{bigq}) and big integers
  (\code{biqz});  consequently, \code{\link{range}()} works as well.

  Similarly, S4 methods are provided for \code{\link{which.min}()} and
  \code{\link{which.max}()}.
}

\arguments{
  \item{x}{a \dQuote{big integer} (\code{bigz}) or \dQuote{big rational}
    (\code{bigq}) vector.}
  \item{...}{numeric arguments}
  \item{na.rm}{a logical indicating whether missing values should be removed.}
}

\value{
  an object of class \code{"\link{bigz}"} or \code{"\link{bigq}"}.
}

\author{Antoine Lucas}% Martin Maechler: which.min & which.max

\examples{
 x <- as.bigz(1:10)
 max(x)
 min(x)
 range(x) # works correctly via default method
 x <- x[c(7:10,6:3,1:2)]
 which.min(x) ## 9
 which.max(x) ## 4

 Q <- as.bigq(1:10, 3)
 max(Q)
 min(Q)
 (Q <- Q[c(6:3, 7:10,1:2)])
 stopifnot(which.min(Q) == which.min(asNumeric(Q)),
           which.max(Q) == which.max(asNumeric(Q)))

stopifnot(range(x) == c(1,10), 3*range(Q) == c(1,10))
}

\seealso{
  \code{\link{max}} etc in \pkg{base}.
}

\keyword{arith}

