\name{summary.wst}
\alias{summary.wst}
\title{Print out some basic information associated with a wst object}
\usage{
\method{summary}{wst}(object, \dots)
}
\arguments{
\item{object}{The object to print a summary about}
\item{\dots}{Other arguments}
}
\description{
Prints out the number of levels, the length of the original vector
from which the object came, the type of wavelet filter associated with
the decomposition, and the date of production.
}
\details{
Description says all
}
\value{
Nothing
}
\seealso{\code{\link{wst}}}
\examples{
vwst <- wst(rnorm(32))
summary(vwst)
#Levels:  5 
#Length of original:  32 
#Filter was:  Daub cmpct on least asymm N=10 
#Date:  Mon Mar  8 21:56:12 2010 
}
\keyword{print}
\author{G P Nason}
