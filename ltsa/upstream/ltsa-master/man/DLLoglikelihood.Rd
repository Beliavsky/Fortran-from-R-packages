\name{DLLoglikelihood}
\alias{DLLoglikelihood}
\title{Durbin-Levinsion Loglikelihood }
\description{
  The Durbin-Levinsion algorithm is used for the computation of the exact
  loglikelihood function.
}
\usage{
DLLoglikelihood(r, z, useC = TRUE)
}

\arguments{
  \item{r}{autocovariance or autocorrelation at lags 0,...,n-1, where n is length(z) }
  \item{z}{time series data}
  \item{useC}{ TRUE, use compiled C, otherwise R }
}
\details{
  The concentrated loglikelihood function may be written Lm(beta) = -(n/2)*log(S/n)-0.5*g,
  where beta is the parameter vector, n is the length of the time series, S=z'M z,
  z is the mean-corrected time series, M is the inverse of the covariance matrix setting
  the innovation variance to one and g=-log(det(M)).
  This method was given in Li (1981) for evaluating the loglikelihood function
  in the case of the fractionally differenced white noise.
}
\value{
The loglikelihood concentrated over the parameter for the innovation
variance is returned.
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
\note{ 
The purpose of this function is to provide a check on the TrenchLoglikelihod function. 
Completely different algorithms are used in each case but the numerical values should
agree.
}
\author{ A.I. McLeod }
\seealso{ 
\code{\link{TrenchLoglikelihood}} 
}
\examples{
#Example 1
#compute loglikelihood for white noise
z<-rnorm(100)
DLLoglikelihood(c(1,rep(0,length(z)-1)), z)

#Example 2
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

#Example 3
\dontrun{
#Compare direct evaluation of AR(1) loglikelihood with DL method
#First define the exact concentrated loglikelihood function for AR(1)
AR1Loglikelihood <-function(phi,z){
n<-length(z)
S<-(z[1]^2)*(1-phi^2) + sum((z[-1]-phi*z[-n])^2)
0.5*log(1-phi^2)-(n/2)*log(S/n)
}
#Next run script to compare numerically the loglikelihoods.
#They should be identical.
phi<-0.8
n<-200
z<-arima.sim(list(ar=phi), n=n)
phis<-seq(0.1,0.95,0.05)
ansAR1<-ansDL<-numeric(length(phis))
for (i in 1:length(phis)) {
    ansAR1[i] <- AR1Loglikelihood(phis[i],z)
    r<-(1/(1-phis[i]^2))*phis[i]^(0:(n-1))
    ansDL[i] <- DLLoglikelihood(r,z,useC=FALSE)
}
ans<-matrix(c(ansDL,ansAR1),ncol=2)
dimnames(ans)<-list(phis, c("DL-method","AR1-method"))
}

#Example 4
\dontrun{
#compare timings. See (McLeod, Yu, Krougly, Table 8).
 n<-5000
 ds<-c(-0.45, -0.25, -0.05, 0.05, 0.25, 0.45)
 tim<-matrix(numeric(3*length(ds)),ncol=3)
 for (i in 1:length(ds)){
    d<-ds[i]
    alpha <- 1-2*d #equivalent hyperbolic autocorrelation
    r <- (1/(1:n))^alpha
    z<-DLSimulate(n,r)
    tim1a<-system.time(LL1<-DLLoglikelihood(r,z))[1]
    tim1b<-system.time(LL1<-DLLoglikelihood(r,z,useC=FALSE))[1]
    tim2<-system.time(LL2<-TrenchLoglikelihood(r,z))[1]
    tim[i,]<-c(tim1a,tim1b, tim2)
    }
 dimnames(tim)<-list(ds, c("DL-C","DL-R","Trench"))
 tim
}
}
\keyword{ts }
