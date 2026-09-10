\name{apply}
\alias{apply}
\alias{apply.bigz}
\alias{apply.bigq}
\alias{apply.default}
\title{Apply Functions Over Matrix Margins (Rows or Columns)}
\usage{
\method{apply}{bigz}(X, MARGIN, FUN, \dots)
\method{apply}{bigq}(X, MARGIN, FUN, \dots)
}
\description{
  These are S3 \code{\link{methods}} for \code{\link{apply}()} which we
  re-export as S3 generic function.
  They \dQuote{overload} the \code{apply()} function for big rationals (\code{"bigq"})
  and big integers (\code{"bigz"}).
}

\arguments{
  \item{X}{a matrix of class bigz or bigq, see e.g.,
    \code{\link{matrix.bigz}}.}
  \item{MARGIN}{1: apply function to rows; 2: apply function to columns}
  \item{FUN}{\code{\link{function}} to be applied}
  \item{\dots}{(optional) extra arguments for \code{FUN()}, as e.g., in
    \code{\link{lapply}}.}
}
\value{
  The \code{bigz} and \code{bigq} methods return a vector of class
  \code{"bigz"} or \code{"bigq"}, respectively.
}
\seealso{
  \code{\link{apply}}; \code{\link{lapply}} is used by our
  \code{apply()} method.
}
\author{Antoine Lucas}

\examples{
 x <- as.bigz(matrix(1:12,3))
 apply(x,1,min)
 apply(x,2,max)

 x <- as.bigq(x ^ 3, d = (x + 3)^2)
 apply(x,1, min)
 apply(x,2, sum)
 ## now use the "..." to pass  na.rm=TRUE :
 x[2,3] <- NA
 apply(x,1, sum)
 apply(x,1, sum, na.rm = TRUE)
}
\keyword{arith}

