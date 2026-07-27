\name{Auto.VR}
\alias{Auto.VR}
\title{ Automatic Variance Ratio Test }
\description{
  A variance ratio test with holding period value chosen by a data dependent procedure
}
\usage{
Auto.VR(y)
}
\arguments{
  \item{y}{ financial return time series }
}

\value{
 
 
  \item{stat }{Automatic variance ratio test statistic}
  \item{sum }{1+ weighted sum of autocorrelation up to the optimal order}
}
\references{ Choi, I. 1999, Testing the random walk hypothesis for real exchange rates
              Journal of Applied Econometrics, 14, 293-308.}
\author{ Jae H. Kim}
\note{ R code translated from Choi's GAUSS code
}
\examples{
data(exrates)
y <- exrates$ca                           
nob <- length(y)
r <- log(y[2:nob])-log(y[1:(nob-1)])      
Auto.VR(r)
}
\keyword{  htest  }
