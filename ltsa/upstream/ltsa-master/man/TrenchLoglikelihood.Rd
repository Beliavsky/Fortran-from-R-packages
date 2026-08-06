\name{TrenchLoglikelihood}
\alias{TrenchLoglikelihood}
\title{Loglikelihood function of stationary time series
  using Trench algorithm}
\description{
  The Trench matrix inversion algorithm is used to compute the
  exact concentrated loglikelihood function.
}
\usage{TrenchLoglikelihood(r, z)}
\arguments{
  \item{r}{autocovariance or autocorrelation at lags 0,...,n-1, where n is length(z) }
  \item{z}{time series data}
}
\details{
  The concentrated loglikelihood function may be written Lm(beta) = -(n/2)*log(S/n)-0.5*g,
  where beta is the parameter vector, n is the length of the time series, S=z'M z,
  z is the mean-corrected time series, M is the inverse of the covariance matrix setting
  the innovation variance to one and g=-log(det(M)).
}
\value{
The loglikelihood concentrated over the parameter for the innovation
variance is returned.
}
\references{ 
McLeod, A.I., Yu, Hao, Krougly, Zinovi L.  (2007).
Algorithms for Linear Time Series Analysis,
Journal of Statistical Software.
}
\author{ A.I. McLeod }
\seealso{ \code{\link{DLLoglikelihood}} }

\examples{
#compute loglikelihood for white noise
z<-rnorm(100)
TrenchLoglikelihood(c(1,rep(0,length(z)-1)), z)


#simulate a time series and compute the concentrated loglikelihood using DLLoglikelihood and
#compare this with the value given by TrenchLoglikelihood.
phi<-0.8
n<-200
r<-phi^(0:(n-1))
z<-arima.sim(model=list(ar=phi), n=n)
LD<-DLLoglikelihood(r,z)
LT<-TrenchLoglikelihood(r,z)
ans<-c(LD,LT)
names(ans)<-c("DLLoglikelihood","TrenchLoglikelihood")
}
\keyword{ts }
