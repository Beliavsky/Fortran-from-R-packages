\name{summary.wst2D}
\alias{summary.wst2D}
\title{Print out some basic information associated with a wst2D object}
\usage{
\method{summary}{wst2D}(object, \dots)
}
\arguments{
\item{object}{The object to print a summary about}
\item{\dots}{Other arguments}
}
\description{
Prints out the number of levels, the dimensions  of the original image
from which the object came, the type of wavelet filter associated with
the decomposition, and the date of production.
}
\details{
Description says all
}
\value{
Nothing
}
\seealso{\code{\link{wst2D}}}
\examples{
m <- matrix(rnorm(32*32), nrow=32)
mwst2D <- wst2D(m)
summary(mwst2D)
#Levels:  5 
#Length of original:  32 x 32 
#Filter was:  Daub cmpct on least asymm N=10 
#Date:  Mon Mar  8 21:57:55 2010 
}
\author{G P Nason}
\keyword{print}
