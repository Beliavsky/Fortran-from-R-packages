\name{Wright}
\alias{Wright}
\title{Wright's Rank and Sign Tests}
\description{
  The function returns R1, R2 and S1 tests statistics detailed in Wright (2000)
}
\usage{
Wright(y, kvec)
}
\arguments{
  \item{y}{ a vector of time series, typically  financial return }
  \item{kvec}{ a vector of holding periods }
}
\details{Nonparametric tests
 
}
\value{
  
  \item{Holding.Period }{holding periods used}
  \item{R1.test }{rank test R1}
  \item{R2.test }{rank test R2}
  \item{S1.test }{sign test S1}
}
\references{ WRIGHT,J.H.,2000,Alternative Variance-Ratio Tests Using Ranks and Signs, Journal of Business & Economic Statistics, 18, 1-9.
}
\author{ Jae H. Kim }


\examples{
data(exrates)
y <- exrates$ca                    
nob <- length(y)
r <- log(y[2:nob])-log(y[1:(nob-1)])  
kvec <- c(2,5,10)
Wright(r,kvec) 
}

\keyword{ htest   }