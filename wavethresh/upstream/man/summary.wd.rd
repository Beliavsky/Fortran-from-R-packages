\name{summary.wd}
\alias{summary.wd}
\title{Print out some basic information associated with a wd object}
\usage{
\method{summary}{wd}(object, \dots)
}
\arguments{
\item{object}{The object to print a summary about}
\item{\dots}{Other arguments}
}
\description{
Prints out the number of levels, the length of the original vector
from which the object came, the type of wavelet filter associated with
the decomposition, the type of boundary handling, the transform type
and the date of production.
}
\details{
Description says all
}
\value{
Nothing
}
\seealso{\code{\link{wd}}}
\examples{
vwd <- wd(1:8)
summary(vwd)
#Levels:  3 
#Length of original:  8 
#Filter was:  Daub cmpct on least asymm N=10 
#Boundary handling:  periodic 
#Transform type:  wavelet 
#Date:  Mon Mar  8 21:30:32 2010 
}
\author{G P Nason}
\keyword{print}
