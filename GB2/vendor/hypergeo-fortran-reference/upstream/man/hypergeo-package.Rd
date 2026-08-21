\name{hypergeo-package}
\alias{hypergeo-package}
\docType{package}
\title{
The hypergeometric function
}
\description{
The hypergeometric function for the whole complex plane
}
\details{
\tabular{ll}{
Package: \tab hypergeo\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-04-16\cr
License: \tab GPL\cr
}
The front end function of the package is \code{hypergeo()}: depending on
the value of the parameters, this executes one or more of many
sub-functions.
}
\author{
Robin K. S. Hankin

Maintainer: <hankin.robin@gmail.com>

}
\references{
  M. Abramowitz and I. A. Stegun 1965. \emph{Handbook of
    mathematical functions}.  New York: Dover
}
\keyword{ package }
\examples{
hypergeo(1.1,2.3,1.9 , 1+6i)

options(showHGcalls = TRUE)    # any non-null value counts as TRUE
hypergeo(4.1, 3.1, 5.1, 1+1i)  # shows trace back
options(showHGcalls = FALSE)   # reset


}
