\name{Chow.Denning}
\alias{Chow.Denning}
\title{ Chow-Denning Multiple Variance Ratio Tests }
\description{
  This function returns Chow-Denning test statistics.

  CD1: test for iid series; CD2: test for uncorrelated series with possible heteroskedasticity.
}
\usage{
Chow.Denning(y, kvec)
}
\arguments{
  \item{y}{ a vector of time series, typically  financial return }
  \item{kvec}{ a vector of holding periods }
}


\value{
  
  \item{Holding.Periods }{holding periods used}
  \item{CD1 }{CD1 statistic}
  \item{CD2 }{CD2 statistic}
  \item{Critical.Values_10_5_1_percent }{10 5 1 percent critical values}
}

\note{

See Chow and Denning (1993) for the details of critical value calculation
}

\references{ Chow,K. V., K. C. DENNING, 1993, A Simple Multiple Variance Ratio Test, Journal of Econometrics, 58, 385-401. }
\author{ Jae H. Kim }


\examples{
data(exrates)
y <- exrates$ca                  
nob <- length(y)
r <- log(y[2:nob])-log(y[1:(nob-1)])      
kvec <- c(2,5,10)
Chow.Denning(r,kvec) 
}

\keyword{  htest  }