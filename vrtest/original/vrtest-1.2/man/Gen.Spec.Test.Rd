\name{Gen.Spec.Test}
\alias{Gen.Spec.Test}
\title{ Generalized spectral Test  }
\description{
  Generalized spectral Test
}
\usage{
Gen.Spec.Test(y,B)
}
\arguments{
  \item{y}{ financial return time series }
  \item{B}{ the number of bootstrap iterations, the default is 300}
}

\value{
 
  \item{Pboot}{wild bootstrap p-value of the test}
}
\references{ 
Escanciano, J.C. and Velasco, C., 2006, Generalized Spectral Tests for the martigale Difference Hypothesis, Journal of Econometrics, 134, p151-185.

Charles, A. Darne, O. Kim, J.H. 2011, Small Sample Proeprties of Alternative Tests for Martingale Difference Hypothesis, Economics Letters, 110(2), 151-154.}
\author{ Jae H. Kim}

\examples{
r <- rnorm(100)        
Gen.Spec.Test(r)
}
\keyword{  htest  }
