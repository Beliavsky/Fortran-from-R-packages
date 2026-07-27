\name{Lo.Mac}
\alias{Lo.Mac}
\title{Lo-MacKinlay variance Ratio Tests }
\description{
  The function returns M1 and M2 statistics of Lo and MacKinlay (1998).

  M1: tests for iid series; M2: for uncorrelated series with possible heteroskedasticity.
}
\usage{
Lo.Mac(y, kvec)
}
\arguments{
  \item{y}{ a vector of time series, typically  financial return }
  \item{kvec}{ a vector of holding periods }
}

\value{
  
  \item{Stats }{ M1 and M2 statistics}
 
}
\references{ LO, A. W., and A. C. MACKINLAY (1988): "Stock Market Prices Do Not Follow Random Walks: Evidence from a Simple Specification Test," The Review of Financial Studies, 1, 41-66. }
\author{ Jae H. Kim}

\examples{
data(exrates)
y <- exrates$ca                 
nob <- length(y)
r <- log(y[2:nob])-log(y[1:(nob-1)])       
kvec <- c(2,5,10)
Lo.Mac(r,kvec) 
}

\keyword{ htest   }