\name{DLSimulate}
\alias{DLSimulate}
\title{ Simulate linear time series }
\description{
The Durbin-Levinsion recursions are used to simulate a stationary
time series given an unit innovation sequence and given autocovariance
function. Requires \deqn{O(n^2)} flops.
}
\usage{
DLSimulate(n, r, useC = TRUE, rand.gen = rnorm, ...)
}

\arguments{
  \item{n}{ length of time series to be generated }
  \item{r}{ autocovariances, lags 0, ...,  }
  \item{useC}{ =TRUE, use C interface. Otherwise direct
   computation. }
  \item{rand.gen}{ random number generator to use}
  \item{\dots}{ optional arguments passed to \code{rand.gen} }
}
\details{
See Hipel and McLeod (1994) or
McLeod, Yu and Krougly (2007).
}
\value{
simulated time series of length n
}
\references{ 
McLeod, A.I., Yu, Hao, Krougly, Zinovi L.  (2007).
Algorithms for Linear Time Series Analysis,
Journal of Statistical Software.
}
\author{ A.I. McLeod }

\seealso{ 
\code{\link{DHSimulate}},
\code{\link{SimGLP}},
\code{\link{arima.sim}} 
}
\examples{
#Simulate hyperbolic decay time series 
#with Hurst coefficient, H=0.9
n<-2000
H<-0.9
alpha<-2*(1-H)  #hyperbolic decay parameter
r<-(1/(1:n))^alpha
z<-DLSimulate(n, r)
plot.ts(z)
#can use HurstK function in FGN library to estimate H

}
\keyword{ ts }
\keyword{ datagen }
