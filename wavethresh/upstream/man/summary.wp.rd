\name{summary.wp}
\alias{summary.wp}
\title{Print out some basic information associated with a wp object}
\usage{
\method{summary}{wp}(object, \dots)
}
\arguments{
\item{object}{The object to print a summary about}
\item{\dots}{Other arguments}
}
\description{
Prints out the number of levels, the length of the original vector
from which the object came, the type of wavelet filter associated with
the decomposition.
}
\details{
Description says all
}
\value{
Nothing
}
\seealso{\code{\link{wp}}}
\examples{
vwp <- wp(rnorm(32))
summary(vwp)
#Levels:  5 
#Length of original:  32 
#Filter was:  Daub cmpct on least asymm N=10 
}
\author{G P Nason}
\keyword{print}
