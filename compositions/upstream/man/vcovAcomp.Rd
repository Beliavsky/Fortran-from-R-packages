\name{vcovAcomp}
\alias{vcovAcomp}
\title{Variance covariance matrix of parameters in compositional regression}
\description{
  The variance covariance tensor structured according of linear models
  with ilr(acomp(...)) responses.
}
\usage{
vcovAcomp(object,...)
          }
\arguments{
  \item{object}{ a statistical model}
  \item{\dots}{further optional parameters for \code{vcov}}
}
\value{
  An array with 4 dimensions. The first 2 are the index dimensions of
  the ilr transform. The later 2 are the index of the parameter.
}

\details{
  The prediction error in compositional linear regression models is a
  complicated object. The function should help to organize it. 
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{vcov}} 
}
\examples{
data(SimulatedAmounts)
model <- lm(ilr(sa.groups)~sa.groups.area)
vcovAcomp(model)[,,1,1]
}
\keyword{multivariate}
