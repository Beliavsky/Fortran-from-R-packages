\name{arima.boot}
\alias{arima.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the Bootstrap Estimates of an ARIMA Model} 
\description{
This function bootstraps time series according to the fitted ARMA(p,d,q) model supplied by the fitted object arima.fit, and estimate the same model using the
arima function. Any bootstrap sample that has problem when fitted with the ARIMA model will be omitted from the final results and all error messages will be suppressed. You can check if there is any fitting problem by running the command geterrmessage().}
\usage{
arima.boot(arima.fit, cond.boot = FALSE, is.normal = TRUE, B = 1000, init, ntrans = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{arima.fit}{a fitted object from the arima function (seasonal components not allowed)} 
  \item{cond.boot}{whether or not 
the bootstrap is conditional on the (p+d) initial values; if
it is set true. If false (default), the stationary bootstrap is used. }
  \item{is.normal}{if true (default), errors are normally distributed, 
otherwise
errors are drawn randomly and with replacement from the residuals of the fitted model.
}
  \item{B}{number of bootstrap replicates (1000, default)}
  \item{init}{initial values for the bootstrap; needed if cond.boot=True 
default values are the initial values of the time series of the fitted 
model.}
  \item{ntrans}{number of transient values for the stationary bootstrap. 
Default=100}
}
\value{
a matrix each row of which consists of the coefficient estimates of 
a bootstrap time-series.}
\author{Kung-Sik Chan}
\examples{
\donttest{
data(hare)
arima.hare=arima(sqrt(hare),order=c(3,0,0))
boot.hare=arima.boot(arima.hare,B=50,init=sqrt(hare)[1:3],ntrans=100)
apply(boot.hare,2,quantile, c(.025,.975))
period.boot=apply(boot.hare,1,function(x){
roots=polyroot(c(1,-x[1:3]))
min1=1.e+9
rootc=NA
for (root in roots) {
if( abs(Im(root))<1e-10) next
if (Mod(root)< min1) {min1=Mod(root); rootc=root}
}
if(is.na(rootc)) period=NA else period=2*pi/abs(Arg(rootc))
period
})
hist(period.boot)
quantile(period.boot,c(0.025,.975))
}
}
\keyword{methods}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
