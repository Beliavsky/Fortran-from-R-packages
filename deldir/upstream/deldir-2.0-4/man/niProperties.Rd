\name{niProperties}
\alias{niProperties}
\docType{data}
\title{
    Northern Ireland properties.
}
\description{
    The locations (in longitude and latitude) of a number of properties
    (land holdings) in Northern Ireland.
}
\usage{data("niProperties")}
\format{
  A data frame with 240 observations on the following 2 variables.
  \describe{
    \item{\code{x}}{A numeric vector of longitudes.}
    \item{\code{y}}{A numeric vector of latitudes.}
  }
}
\source{
   These data were kindly provided by Adam Dadvar of the
   \emph{Cartesian Limited} consulting service.\cr
    URL: \code{http://www.cartesian.com}.
}
\examples{
# data(niProperties)
# It is unnecessary to use \code{data} since \code{niProperties} is
# a "first class object".  It is "lazily loaded".
plot(niProperties)
}
\keyword{datasets}
