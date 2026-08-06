\name{DLResiduals}
\alias{DLResiduals}
\title{ Prediction residuals }
\description{
The Durbin-Levison algorithm is used to compute the one-step
prediction residuals.
}
\usage{
DLResiduals(r, z, useC = TRUE, StandardizedQ=TRUE)
}
\arguments{
  \item{r}{ vector of length n containing the autocovariances or autocorrelations at lags 0,...,n-1 }
  \item{z}{ vector of length n, mean-corrected time series data }
  \item{useC}{ if TRUE, the compiled C code is used, otherwise
  the computations are done entirely in R and much slower}
  \item{StandardizedQ}{TRUE, the residuals are divided by their standard deviation or FALSE, the
  raw prediction residuals are computed}
}
\details{
If the model is correct the standardized prediction residuals are approximately NID(0,1) and 
are asymptotically
equivalent to the usual innovation residuals divided by the residual sd. This means that
the usual diagnotic checks, such as the Ljung-Box test may be used.
  }
\value{
Vector of length n containing the residuals
}
\references{ 
W.K. Li (1981). 
Topics in Time Series Analysis. 
Ph.D. Thesis, 
University of Western Ontario.

McLeod, A.I., Yu, Hao, Krougly, Zinovi L.  (2007).
Algorithms for Linear Time Series Analysis,
Journal of Statistical Software.
}
\author{ A.I. McLeod }
\seealso{ \code{\link{DLLoglikelihood}} }
\examples{
# For the AR(1) the prediction residuals and innovation residuals are the same (except for
# t=1).  In this example we demonstrate the equality of these two types of residuals.
#
phi<-0.8
sde<-30
n<-30
z<-arima.sim(n=30,list(ar=phi),sd=sde)
r<-phi^(0:(n-1))/(1-phi^2)*sde^2
e<-DLResiduals(r,z)
a<-numeric(n)
for (i in 2:n)
    a[i]=z[i]-phi*z[i-1]
a<-a/sde
ERR<-sum(abs(e[-1]-a[-1]))
ERR
#
#Simulate AR(1) and compute the MLE for the innovation variance
phi <- 0.5
n <- 2000
sigsq <- 9
z<-arima.sim(model=list(ar=phi), n=n, sd=sqrt(sigsq))
g0 <- sigsq/(1-phi^2)
r <- g0*phi^(0:(n-1))
#comparison of estimate with actual
e<-DLResiduals(r,z,useC=FALSE, StandardizedQ=FALSE)
sigsqHat <- var(e)
ans<-c(sigsqHat,sigsq)
names(ans)<-c("estimate","theoretical")
ans
}
\keyword{ ts }
