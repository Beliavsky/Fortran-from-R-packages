\name{PredictionVariance}
\alias{PredictionVariance}
\title{ Prediction variance}
\description{
The prediction variance of the forecast for lead times l=1,...,maxLead
is computed given theoretical autocovariances. 
}
\usage{
PredictionVariance(r, maxLead = 1, DLQ = TRUE)
}

\arguments{
  \item{r}{ the autocovariances at lags 0, 1, 2, ... }
  \item{maxLead}{ maximum lead time of forecast }
  \item{DLQ}{ Using Durbin-Levinson if TRUE. Otherwise Trench algorithm used. }
}

\details{
Two algorithms are available which 
are described in detail in McLeod, Yu and Krougly (2007).
The default method, DLQ=TRUE, uses the autocovariances provided in r to
determine the optimal linear mean-square error predictor of order
length(r)-1.  
The mean-square error of this predictor is the lead-one error variance.
The moving-average expansion of this model is used to compute any
remaining variances (McLeod, Yu and Krougly, 2007).
With the other Trench algorithm, when DLQ=FALSE, a direct matrix representation
of the forecast variances is used (McLeod, Yu and Krougly, 2007).
The Trench method is exact.  Provided the length of r is large enough,
the two methods will agree.
}

\value{
vector of length maxLead containing the variances
}

\references{ 
McLeod, A.I., Yu, Hao, Krougly, Zinovi L.  (2007).
Algorithms for Linear Time Series Analysis,
Journal of Statistical Software.
}

\author{ A.I. McLeod }

\seealso{ 
\code{\link{predict.Arima}}, 
\code{\link{TrenchForecast}},
\code{\link{exactLoglikelihood}}
}

\examples{
#Example 1. Compare using DL method or Trench method
va<-PredictionVariance(0.9^(0:10), maxLead=10)
vb<-PredictionVariance(0.9^(0:10), maxLead=10, DLQ=FALSE)
cbind(va,vb)
# 
#Example 2. Compare with predict.Arima
#general script, just change z, p, q, ML
z<-sqrt(sunspot.year)
n<-length(z)
p<-9
q<-0
ML<-10
#for different data/model just reset above
out<-arima(z, order=c(p,0,q))
sda<-as.vector(predict(out, n.ahead=ML)$se)
#
phi<-theta<-numeric(0)
if (p>0) phi<-coef(out)[1:p]
if (q>0) theta<-coef(out)[(p+1):(p+q)]
zm<-coef(out)[p+q+1]
sigma2<-out$sigma2
r<-sigma2*tacvfARMA(phi, theta, maxLag=n+ML-1)
sdb<-sqrt(PredictionVariance(r, maxLead=ML))
cbind(sda,sdb)
#
# 
#Example 3. DL and Trench method can give different results
#  when the acvf is slowly decaying. Trench is always
#  exact based on a finite-sample.
L<-5
r<-1/sqrt(1:(L+1))
va<-PredictionVariance(r, maxLead=L)
vb<-PredictionVariance(r, maxLead=L, DLQ=FALSE)
cbind(va,vb) #results are slightly different
r<-1/sqrt(1:(1000)) #larger number of autocovariances
va<-PredictionVariance(r, maxLead=L)
vb<-PredictionVariance(r, maxLead=L, DLQ=FALSE)
cbind(va,vb) #results now agree
}

\keyword{ ts }

