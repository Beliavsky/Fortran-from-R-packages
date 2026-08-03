\name{skewness}
\alias{skewness}
\title{ Skewness of the sample }
\description{
  This function computes skewness of given data.
}
\usage{
skewness(x, na.rm = FALSE)
}
\arguments{
  \item{x}{ a numeric vector, matrix or data frame. }
  \item{na.rm}{ logical. Should missing values be removed? }
}
\author{ Lukasz Komsta }
\seealso{ \code{\link{agostino.test}} }
\examples{
set.seed(1234)
skewness(rnorm(1000))
}
\keyword{univar}

