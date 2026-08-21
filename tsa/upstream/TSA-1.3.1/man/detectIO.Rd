\name{detectIO}
\alias{detectIO}
\title{Innovative Outlier Detection}
\description{
This function serves to detect whether there are any innovative
outliers (IO). It implements the
test statistic \eqn{lambda_{2,t}} proposed by Chang, Chen and Tiao (1988).
}
\usage{
detectIO(object, alpha = 0.05, robust = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a fitted ARIMA model}
  \item{alpha}{family significance level (5\% is the default)
Bonferroni rule is used to control the family error rate.}
  \item{robust}{if true, the noise standard deviation is estimated by
mean absolute residuals times sqrt(pi/2). Otherwise,
it is the estimated by sqrt(sigma2) from the arima fit.}
}
\value{
A list containing the following components:
\item{ind}{the time indices of potential AO}
\item{lambda1}{the corresponding test statistics}
}
\references{ Chang, I.H., Tiao, G.C. and C. Chen (1988). Estimation of Time Series Parameters in the Presence of Outliers. 
Technometrics, 30, 193-204.}
\author{Kung-Sik Chan}
\seealso{ \code{\link{detectIO}}}
\examples{
set.seed(12345)
y=arima.sim(model=list(ar=.8,ma=.5),n.start=158,n=100)
y[10]
y[10]=10
y=ts(y,freq=1,start=1)
plot(y,type='o')
acf(y)
pacf(y)
eacf(y)
m1=arima(y,order=c(1,0,0))
m1
detectAO(m1)
detectAO(m1, robust=FALSE)
detectIO(m1)
}
\keyword{methods}
