\name{harmonic}
\alias{harmonic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Construct harmonic functions for fitting harmonic trend model}
\description{
The function creates a matrix of the first m pairs of harmonic functions for
fitting a harmonic trend 
(cosine-sine trend, Fourier regresssion) models with the response
being x, a time series.  
}
\usage{
harmonic(x, m = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a time series}
  \item{m}{the number of pairs of harmonic functions to be created; 2m must be
less than or equal to the frequency of x}
}
\value{a matrix consisting of \eqn{\cos(2k \pi t), \sin(2k \pi t), k=1,2,...,m,}
 excluding any zero functions.
}
\author{Kung-Sik Chan}
\seealso{ \code{\link{season}}}
\examples{
data(tempdub)
# first creates the first pair of harmonic functions and then fit the model
har.=harmonic(tempdub,1)
model4=lm(tempdub~har.)
summary(model4)
}
\keyword{methods}
