\name{summary.wd3D}
\alias{summary.wd3D}
\title{Print out some basic information associated with a wd3D object}
\usage{
\method{summary}{wd3D}(object, \dots)
}
\arguments{
\item{object}{The object to print a summary about}
\item{\dots}{Other arguments}
}
\description{
Prints out the number of levels, the type of wavelet filter associated with
the decomposition, and the date of production.
}
\details{
Description says all
}
\value{
Nothing
}
\seealso{\code{\link{wd3D}}}
\examples{
test.data.3D <- array(rnorm(8*8*8), dim=c(8,8,8))
tdwd3D <- wd3D(test.data.3D)
summary(tdwd3D)
#Levels:  3 
#Filter number was:  10 
#Filter family was:  DaubLeAsymm 
#Date:  Mon Mar  8 21:48:00 2010 
}
\author{G P Nason}
\keyword{print}
