\name{Adjust.thin}
\alias{Adjust.thin}

\title{ Adjustment for thinly-traded returns }
\description{
  The adjustment based on AR(1) fitting as proposed by Miller et al. (1994)
}
\usage{
Adjust.thin(y)
}
\arguments{
  \item{y}{ financial return time series }
}

\value{
  
  Adjusted return
}

\references{ 
Miller et al. (1994), Mean Reversion of Standard & Poor's 500 Index Base Changes: Arbitrage Induced or Statistical Illusion
  Journal of Finance, XLIX, 479-513.}
\author{ Jae H. Kim }


\examples{
data(exrates)
y <- exrates$ca                           
nob <- length(y)
r <- log(y[2:nob])-log(y[1:(nob-1)]) 
Adjust.thin(r)
}
\keyword{ htest }
