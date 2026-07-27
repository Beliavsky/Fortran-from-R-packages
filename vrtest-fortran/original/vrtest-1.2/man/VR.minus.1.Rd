\name{VR.minus.1}
\alias{VR.minus.1}
\title{ Absolute Value of (VR - 1) }
\description{
  This value is sometimes used to measure the degree of market efficiency
}
\usage{
VR.minus.1(y, kvec)
}
\arguments{
  \item{y}{ financial return time series }
  \item{kvec}{ a vector of holding periods }
}

\value{
  
  \item{VR.auto }{ the value of VR-1 with automatic selection of holding vectors}
  \item{Holding.Peiods }{the vector of holding periods}
  \item{VR.kvec }{the values of VR-1 for the chosen holding periods}
}

\author{ Jae H. Kim  }
\note{ see Auto.VR function for automatic selection of holding periods
}

\examples{
data(exrates)
y <- exrates$ca             
nob <- length(y)
kvec <- c(2,5,10)
r <- log(y[2:nob])-log(y[1:(nob-1)])          
VR.minus.1(r,kvec)
}
\keyword{htest }

