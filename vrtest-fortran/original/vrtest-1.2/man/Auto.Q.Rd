\name{Auto.Q}
\alias{Auto.Q}

\title{ Automatic Portmanteau Test }
\description{
  A robustified portmanteau test with automatic lag selection 
}
\usage{
Auto.Q(y,lags)
}
\arguments{
  \item{y}{ financial return time series }
  \item{lags}{ maximum lag value, the default is 10}
}

\value{
 
 
  \item{Stat }{Automatic portmanteau test statistic}
  \item{Pvalue}{p-value of the test}
}
\references{ 
Escanciano, J.C., Lobato, I.N. 2009a. An automatic portmanteau test for serial correlation. Journal of Econometrics 151, 140-149.

Charles, A. Darne, O. Kim, J.H. 2011, Small Sample Proeprties of Alternative Tests for Martingale Difference Hypothesis, Economics Letters, 110(2), 151-154.
}
\author{ Jae H. Kim}

\examples{
data(exrates)
y <- exrates$ca                          
nob <- length(y)
r <- log(y[2:nob])-log(y[1:(nob-1)])    
Auto.Q(r)
}
\keyword{  htest  }
