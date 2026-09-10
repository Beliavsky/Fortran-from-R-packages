\name{formatN}
\title{Format Numbers Keeping Classes Distinguishable}
\alias{formatN}
\alias{formatN.default}
\alias{formatN.integer}
\alias{formatN.double}
\alias{formatN.bigz}
\alias{formatN.bigq}
\description{
  Format (generalized) numbers in a way that their \code{\link{class}}es
  are distinguishable.  Contrary to \code{\link{format}()} which uses a
  common format for all elements of \code{x}, here, each entry is
  formatted individually.
}
\usage{
formatN(x, ...)
\S3method{formatN}{default}(x, ...)
\S3method{formatN}{integer}(x, ...)
\S3method{formatN}{double}(x, ...)
\S3method{formatN}{bigz}(x, ...)
\S3method{formatN}{bigq}(x, ...)
}
\arguments{
  \item{x}{any \R object, typically \dQuote{number-like}.}
  \item{\dots}{potentially further arguments passed to methods.}
}
\value{
  a character vector of the same \code{\link{length}} as \code{x},
  each entry a representation of the corresponding entry in \code{x}.
}
\author{Martin Maechler}
\seealso{
  \code{\link{format}}, including its (sophisticated) default method;
  \code{\link{as.character}}.

}
\examples{
## Note that each class is uniquely recognizable from its output:
formatN(    -2:5)# integer
formatN(0 + -2:5)# double precision
formatN(as.bigz(-2:5))
formatN(as.bigq(-2:5, 4))
}
\keyword{character}
\keyword{print}
