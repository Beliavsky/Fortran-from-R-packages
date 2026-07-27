\name{Spec.shape}
\alias{Spec.shape}
\title{ Spectral shape tests for random walk}
\description{
  Spectral Shape tests proposed by Durlauf (1991) and Choi (1999)
}
\usage{
Spec.shape(x)
}
\arguments{
  \item{x}{ financial return time series  }
}

\value{
  \item{AD }{ Anderson-Darling statistic}
  \item{CVM }{ Cramer-von Mises statistic}
  \item{M }{ Mellows statistic}
}
\references{ 
Choi, I. 1999, Testing the random walk hypothesis for real exchange rates, Journal of Applied Econometrics, 14, 293-308.
Durlauf, S. N., 1991, Spectral based testing of the martingale hypothesis, Journal of Econometrics, 50, 355-376.
 }
\author{ Jae H. Kim }
\note{  Traslated from Choi's Gauss codes
}

\examples{
data(exrates)
y <- exrates$ca                   
nob <- length(y)
r <- log(y[2:nob])-log(y[1:(nob-1)])     
Spec.shape(r)
}

\keyword{ htest }

