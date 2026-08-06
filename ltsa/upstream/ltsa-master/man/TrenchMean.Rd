\name{TrenchMean}
\alias{TrenchMean}
\title{ Exact MLE for mean given the autocorrelation function }
\description{
  Sometimes this is also referred to as the BLUE.
}
\usage{
TrenchMean(r, z)
}
\arguments{
  \item{r}{ vector of autocorrelations or autocovariances of length n }
  \item{z}{ time series data vector of length n }
}

\value{
the estimate of the mean
}

\references{ 
McLeod, A.I., Yu, Hao, Krougly, Zinovi L.  (2007).
Algorithms for Linear Time Series Analysis,
Journal of Statistical Software.
}

\author{ A.I. McLeod }

\note{ An error is given if r is not a postive-definite sequence or
if the lengths of \code{r} and \code{z} are not equal.
}

\seealso{ \code{\link{TrenchInverse}} }
\examples{
#compare BLUE and sample mean
phi<- -0.9
a<-rnorm(100)
z<-numeric(length(a))
phi<- -0.9
n<-100
a<-rnorm(n)
z<-numeric(n)
mu<-100
sig<-10
z[1]<-a[1]*sig/sqrt(1-phi^2)
for (i in 2:n)
	z[i]<-phi*z[i-1]+a[i]*sig
z<-z+mu
r<-phi^(0:(n-1))
meanMLE<-TrenchMean(r,z)
meanBLUE<-mean(z)
ans<-c(meanMLE, meanBLUE)
names(ans)<-c("BLUE", "MLE")
ans
}
\keyword{ ts }
