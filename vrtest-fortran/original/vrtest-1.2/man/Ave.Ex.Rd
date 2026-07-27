\name{Ave.Ex}
\alias{Ave.Ex}
\title{ Average Exponential Tests}
\description{
  Average exponential tests of Andrews and Ploberger (1996)
}
\usage{
Ave.Ex(y)
}
\arguments{
  \item{y}{ financial return time series}
}

\value{
  
  \item{ Ex.LM }{LM test}
  \item{ Ex.LR }{LR test}
}
\references{ Choi, I. 1999, Testing the random walk hypothesis for real exchange rates, Journal of Applied Econometrics, 14, 293-308. }

\author{ Jae H. Kim }
\note{ Traslated from Choi's Gauss codes
}

\examples{
data(exrates)
y <- exrates$ca                      
nob <- length(y)
r <- log(y[2:nob])-log(y[1:(nob-1)])   
Ave.Ex(r)
}
\keyword{ htest }

