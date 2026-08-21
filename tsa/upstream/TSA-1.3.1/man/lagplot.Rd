\name{lagplot}
\alias{lagplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Lagged Regression Plot}
\description{
Computes and plots the nonparametric regression function of a time series against
its various lags. 
}
\usage{
lagplot(x, lag.max = 6, deg = 1, nn = 0.7, method = c("locfit", "gam", "both")[1])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{time series}
  \item{lag.max}{maximum lag}
  \item{deg}{degree of local polynomial, needed only for the locfit method}
  \item{nn}{fraction of nearest data contained in a window, needed only for the
locfit method}
  \item{method}{Two methods for nonparametric estimation: "locfit" is
the default which 
uses the local polynomial approach via the locfit library to estimate the
conditional mean function of \eqn{E(X_t|X_{t-k}=x)} for \eqn{1\le k \le lag.max}; 
Another method is GAM, via the mgcv library.} 
}
\value{
Side effects: The nonparametric lagged regression functions are plotted lag by
lag, with the raw data superimposed on the plots. 
}
\references{"Time Series Analysis, with Applications in R" by J.D. Cryer and K.S. Chan}
\author{Kung-Sik Chan}
\examples{
set.seed(2534567)
par(mfrow=c(3,2))
y=arima.sim(n=61,model=list(ar=c(1.6,-0.94),ma=-0.64))
# lagplot(y)
}
\keyword{methods}
