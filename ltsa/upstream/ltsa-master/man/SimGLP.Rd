\name{SimGLP}
\alias{SimGLP}

\title{ Simulate GLP given innovations}
\description{
Simulates a General Linear Time Series that can have nonGaussian innovations.
It uses the FFT so it is O(N log(N)) flops where N=length(a) and N is assumed
to be a power of 2. 
The R function \code{convolve} is used which implements the FFT.

}
\usage{
SimGLP(psi, a)
}
\arguments{
  \item{psi}{ vector, length Q, of MA coefficients starting with 1. }
  \item{a}{ vector, length Q+n, of innovations, where n is the length of time series
  to be generated. }
}
\details{

  \deqn{ z_t = \sum_{k=0}^Q psi_k a_{t-k} }
  where \eqn{t=1,\ldots,n} and the innovations
  $a_t, t=1-Q, \ldots, 0, 1, \ldots, n$ are
  given in the input vector a.
  
  Since \code{convolve} uses the FFT this is faster than direct computation.
}
\value{
vector of length n, where n=length(a)-length(psi)
}

\author{ A.I. McLeod }

\seealso{ 
\code{\link{convolve}}, 
\code{\link{arima.sim}} 
 }

\examples{
#Simulate an AR(1) process with parameter phi=0.8 of length n=100 with
#  innovations from a t-distribution with 5 df and plot it.
#
phi<-0.8
psi<-phi^(0:127)
n<-100
Q<-length(psi)-1
a<-rt(n+Q,5)
z<-SimGLP(psi,a)
z<-ts(z)
plot(z)
}
\keyword{ ts }
\keyword{ datagen }
