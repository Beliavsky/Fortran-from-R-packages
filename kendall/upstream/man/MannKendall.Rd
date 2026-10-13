\name{MannKendall}
\alias{MannKendall}
\title{ Mann-Kendall trend test }
\description{
This is a test for monotonic trend in a time series z[t] based on
the Kendall rank correlation of z[t] and t.
}

\details{
The test was suggested by Mann (1945) and has been extensively used with
environmental time series (Hipel and McLeod, 2005).
For autocorrelated time series, the block bootstrap may be used
to obtain an improved signficance test.
}

\usage{
MannKendall(x)
}
\arguments{
  \item{x}{ a vector of data, often a time series }
}
\value{
A list with class Kendall.
  \item{tau}{ Kendall's tau statistic}
  \item{sl }{two-sided p-value}
  \item{S}{Kendall Score}
  \item{D }{denominator, tau=S/D}
  \item{varS }{variance of S}
Generic function print.Kendall and summary.Kendall are provided to
print the output.
}

\note{ If you want to use the output from MannKendall,
save the result as in res<-MannKendall(x,y) and then select from the list res the
value(s) needed. }

\references{
Davison, A.C. and Hinkley, D.V. (1997) Bootstrap Methods and Their Application. Cambridge University Press.

Hipel, K.W. and McLeod, A.I., (2005).
Time Series Modelling of Water Resources and Environmental Systems.
Electronic reprint of our book orginally published in 1994.
\url{http://www.stats.uwo.ca/faculty/aim/1994Book/}.

Mann, H.B. (1945), Nonparametric tests against trend, Econometrica, 13, 245-259.
}

\author{ A.I. McLeod, aimcleod@uwo.ca }


\seealso{ \code{\link{Kendall}}, \code{\link{SeasonalMannKendall}} \code{\link[boot]{tsboot}}  }

\examples{
# Annual precipitation entire Great Lakes
# The time series plot with lowess smooth suggests an upward trend
# The autocorrelation in this data does not appear significant.
# The Mann-Kendall trend test confirms the upward trend.
data(PrecipGL)
plot(PrecipGL)
lines(lowess(time(PrecipGL),PrecipGL),lwd=3, col=2)
acf(PrecipGL)
MannKendall(PrecipGL)
#
#Use block bootstrap 
library(boot)
data(PrecipGL)
MKtau<-function(z) MannKendall(z)$tau
tsboot(PrecipGL, MKtau, R=500, l=5, sim="fixed")
#
# Deseasonalize a monthly time series and use the block bootstrap
library(boot)
data(manaus)
z<-matrix(manaus, ncol=12, byrow=12)
zm<-apply(z, MARGIN=2, FUN=mean)
zs<-apply(z, MARGIN=2, FUN=sd)
z2<-sweep(z, MARGIN=2, STATS=zm) #subtract monthly means
z3<-sweep(z2, MARGIN=2, STATS=zs, FUN="/") #divide by monthly sd
zds<-c(t(z3)) 
attributes(zds)<-attributes(manaus)
plot(zds)
#do Mann-Kendall trend test
MannKendall(zds)
#check robustness by applying block bootstrap
MKtau<-function(z) MannKendall(z)$tau
tsboot(zds, MKtau, R=500, l=12, sim="fixed")
#the significance level has dropped!
#it is significant now at about 2 sd
}


\keyword{ nonparametric }
\keyword{ ts }
