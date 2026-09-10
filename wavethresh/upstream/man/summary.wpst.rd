\name{summary.wpst}
\alias{summary.wpst}
\title{Print out some basic information associated with a wpst object}
\usage{
\method{summary}{wpst}(object, \dots)
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
\seealso{\code{\link{wpst}}}
\examples{
vwpst <- wpst(rnorm(32))
summary(vwpst)
#Levels:  5 
#Length of original:  32 
#Filter was:  Daub cmpct on least asymm N=10 
#Date:  Mon Mar  8 21:54:47 2010 
}
\author{G P Nason}
\keyword{print}
