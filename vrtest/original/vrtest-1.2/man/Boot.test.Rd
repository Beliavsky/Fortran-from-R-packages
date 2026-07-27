\name{Boot.test}
\alias{Boot.test}
\title{ Bootstrap Variance Ratio Tests }
\description{
  This function returns bootstrap p-values of the Lo-MacKilay (1988) and Chow-Denning (1993) tests.

  Users can choose between iid bootstrap and wild bootstrap
}
\usage{
Boot.test(y, kvec, nboot, wild, prob=c(0.025,0.975))
}
\arguments{
  \item{y}{ a vector of time series, typically  financial return}
  \item{kvec}{ a vector of holding periods }
  \item{nboot}{ the number of bootstrap iterations }
  \item{wild}{ "No" for iid bootstrap, "Normal" for the wild bootstrap using the standard normal distribution, "Mammen" for the wild bootstrap using Mammen's two point distribution,
               "Rademacher" for the wild bootstrap using Rademacher's two point distribution}
   \item{prob}{probability limits for confidence intervals}
}

\value{
  
  \item{Holding.Period }{holding periods used}
  \item{LM.pval }{Bootstrap p-values for the Lo-MacKinlay tests}
  \item{CD.pval }{Bootstrap p-value for the Chow-Denning test}
  \item{CI }{Confidence Intervals for Lo-Mackinlay tests from Bootstrap distribution}
}
\references{ Kim, J.H., 2006, Wild Bootstrapping Variance Ratio Tests. Economics Letters, 92, 38-43. }
\author{ Jae H. Kim }


\examples{
data(exrates)
y <- exrates$ca                   
nob <- length(y)
r <- log(y[2:nob])-log(y[1:(nob-1)])    
kvec <- c(2,5,10)
Boot.test(r,kvec,nboot=500,wild="Normal")
}
\keyword{  htest  }
