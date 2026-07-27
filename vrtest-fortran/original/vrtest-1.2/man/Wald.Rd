\name{Wald}
\alias{Wald}
\title{ Wald Test of Richardson and Smith (1991) }
\description{
  This function returns the Wald test statistic with critical values
}
\usage{
Wald(y, kvec)
}
\arguments{
  \item{y}{ a vector of time series, typically  financial return }
  \item{kvec}{ a vector of holding periods }

  
}
\value{
  
  \item{Holding.Periods}{holding periods used}
  \item{Wald.stat}{Wald test statistic}
  \item{Critical.Values_10_5_1_percent}{10 5 and 1 percent critical values}
}
\references{ Richardson, M., T. Smith, 1991, "Tests of Financial Models in the Presence of Overlapping Observations," The Review Financial Studies, 4, 227-254. }
\author{ Jae H. Kim }
\note{ The statistic asymptotically follows the chi-squared distribution with the degrees of freedom same as the number of holding periods used}
\examples{
data(exrates)
y <- exrates$ca                        
nob <- length(y)
r <- log(y[2:nob])-log(y[1:(nob-1)])     
kvec <- c(2,5,10)
Wald(r,kvec) 
}
\keyword{  htest  }