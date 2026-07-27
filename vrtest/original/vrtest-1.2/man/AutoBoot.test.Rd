\name{AutoBoot.test}
\alias{AutoBoot.test}
\title{ Wild Bootstrapping of Automatic Variance Ratio Test }
\description{
  This function returns wild bootstrap test results for the Automatic Variance Ratio Test of Choi (1999)
}
\usage{
AutoBoot.test(y, nboot, wild,prob=c(0.025,0.975))
}
\arguments{
  \item{y}{ a vector of time series, typically  financial return}
  \item{nboot}{ the number of bootstrap iterations }
  \item{wild}{ "Normal" for the wild bootstrap using the standard normal distribution, "Mammen" for the wild bootstrap using Mammen's two point distribution,
               "Rademacher" for the wild bootstrap using Rademacher's two point distribution}
  \item{prob}{probability limits for confidence intervals}
}

\value{
  
  \item{test.stat }{Automatic variance ratio test statistic}
  \item{VRsum }{1+ weighted sum of autocorrelation up to the optimal order}
  \item{pval }{Wild Bootstrap p-value for the test}
  \item{CI }{Confidence Intervals for the test statistic from Bootstrap distribution}
  \item{CI.VRsum }{Confidence Intervals for the VRsum from Bootstrap distribution}
}
\references{

Kim, J. H., 2009, Automatic Variance Ratio Test under Conditional Heteroskedascity, Finance Research Letters, 6(3), 179-185. 

Charles, A. Darne, O. Kim, J.H. 2011, Small Sample Proeprties of Alternative Tests for Martingale Difference Hypothesis, Economics Letters, 110(2), 151-154.}
\author{ Jae H. Kim }


\examples{                                
r <- rnorm(100)          
AutoBoot.test(r,nboot=500,wild="Normal")
}
\keyword{  htest  }

