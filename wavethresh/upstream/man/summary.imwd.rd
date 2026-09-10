\name{summary.imwd}
\alias{summary.imwd}
\title{Print out some basic information associated with an imwd object}
\usage{
\method{summary}{imwd}(object, \dots)
}
\arguments{
\item{object}{The object to print a summary about}
\item{\dots}{Other arguments}
}
\description{
Prints out the number of levels, the dimensions of the original image
from which the object came, the type of wavelet filter associated with
the decomposition, the type of boundary handling.
}
\details{
Description says all
}
\value{
Nothing
}
\seealso{\code{\link{imwd}}, \code{\link{threshold.imwd}}}
\examples{
m <- matrix(rnorm(32*32),nrow=32)
mimwd <- imwd(m)
summary(mimwd)
#UNcompressed image wavelet decomposition structure
#Levels:  5 
#Original image was 32 x 32  pixels.
#Filter was:  Daub cmpct on least asymm N=10 
#Boundary handling:  periodic 
}
\author{G P Nason}
\keyword{print}
