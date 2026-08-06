\name{TrenchForecast}
\alias{TrenchForecast}
\title{ Minimum Mean Square Forecast}
\description{
Given time series of length n+m, the forecasts for lead times k=1,...,L are 
computed starting with forecast origin at time t=n and continuing up to t=n+m.
The input time series is of length n+m.
For purely out-of-sample forecasts we may take n=length(z).
Note that the parameter m is inferred using the fact that m=length(z)-n.
}
\usage{
TrenchForecast(z, r, zm, n, maxLead, UpdateAlgorithmQ = TRUE)
}

\arguments{
  \item{z}{time series data, length n+m  }
  \item{r}{autocovariances of length(z)+L-1 or until damped out  }
  \item{zm}{mean parameter in model  }
  \item{n}{forecast origin, n }
  \item{maxLead}{ =L, the maximum lead time }
  \item{UpdateAlgorithmQ}{ = TRUE, use efficient update method, otherwise if
  UpdateAlgorithmQ=FALSE, the direct inverse matrix is computed each time}
}

\details{
The minimum mean-square error forecast of z[N+k] given time series data z[1],...,z[N]
is denoted by \eqn{z_N(k)}{z_N(k)}, where N is called the forecast origin and k is
the lead time.
This algorithm computes a table for 
\eqn{z_N(k), N=n,\dots,n+m; k=1,\ldots,m}{z_N(k),N=n,...,n+m; k=1,...,m} 
The minimum mean-square error forecast is simply the conditional expectation
of \eqn{z_{N+k}}{z_{N+k}} given the time series up to including time \eqn{t=N}{t=N}.
This conditional expectation works out to the same thing as the conditional
expectation in an appropriate multivariate normal distribution -- even
if no normality assumption is made.
See McLeod, Yu, Krougly (2007, eqn. 8).
Similar remarks hold for the variance of the forecast.
An error message is given if length(r) < n + L -1.
}

\value{
A list with components
  \item{Forecasts }{matrix with m+1 rows and maxLead columns with the forecasts}
  \item{SDForecasts }{matrix with m+1 rows and maxLead columns with the sd of the forecasts}
}
\references{
McLeod, A.I., Yu, Hao, Krougly, Zinovi L.  (2007).
Algorithms for Linear Time Series Analysis,
Journal of Statistical Software.
}
\author{ A.I. McLeod }

\note{
An error message is given if r is not a pd sequence, that is, the Toeplitz
matrix of r must be pd.
This could occur if you were to approximate a GLP which is near the stationary
boundary by a MA(Q) with Q not large enough.
In the bootstrap simulation experiment reported in our paper McLeod, Yu and Krougly (2007)
we initially approximated the FGN autocorrelations by setting them to zero after lag
553 but in this case the ARMA(2,1) forecasts were always better.
When we used all required lags of the acvf then the FGN forecasts were better as
we expected. 
From this experience, we don't recommend setting high-order acf lags to zero unless
the values are in fact very small.
}

\seealso{ \code{\link{TrenchInverse}} }

\examples{
#Example 1. Compare TrenchForecast and predict.Arima
#general script, just change z, p, q, ML
z<-sqrt(sunspot.year)
n<-length(z)
p<-9
q<-0
ML<-10
#for different data/model just reset above
out<-arima(z, order=c(p,0,q))
Fp<-predict(out, n.ahead=ML)

phi<-theta<-numeric(0)
if (p>0) phi<-coef(out)[1:p]
if (q>0) theta<-coef(out)[(p+1):(p+q)]
zm<-coef(out)[p+q+1]
sigma2<-out$sigma2
#r<-var(z)*ARMAacf(ar=phi, ma=theta, lag.max=n+ML-1)
#When r is computed as above, it is not identical to below
r<-sigma2*tacvfARMA(phi, theta, maxLag=n+ML-1)
F<-TrenchForecast(z, r, zm, n, maxLead=ML)
#the forecasts are identical using tacvfARMA
#    
#Example 2. Compare AR(1) Forecasts.  Show how 
#Forecasts from AR(1) are easily calculated directly. 
#We compare AR(1) forecasts and their sd's.
#Define a function for the AR(1) case
AR1Forecast <- function(z,phi,n,maxLead){
        nz<-length(z)
        m<-nz-n
        zf<-vf<-matrix(numeric(maxLead*m),ncol=maxLead)
        zorigin<-z[n:nz]
        zf<-outer(zorigin,phi^(1:maxLead))
        vf<-matrix(rep(1-phi^(2*(1:maxLead)),m+1),byrow=TRUE,ncol=maxLead)/(1-phi^2)
        list(zf=zf,sdf=sqrt(vf))
        }
#generate AR(1) series and compare the forecasts
phi<-0.9
n<-200
m<-5
N<-n+m
z<-arima.sim(list(ar=phi), n=N)
maxLead<-3
nr<-N+maxLead-1
r<-(1/(1-phi^2))*phi^(0:nr) 
ansT1<-TrenchForecast(z,r,0,n,maxLead)
ansT2<-TrenchForecast(z,r,0,n,maxLead,UpdateAlgorithmQ=FALSE)
ansAR1<-AR1Forecast(z,phi,n,maxLead)
}
\keyword{ ts }
