\name{DL.test}
\alias{DL.test}
\title{ Dominguez-Lobato Test for Martingale Difference Hypothesis }
\description{
  Dominguez-Lobato Test
}
\usage{
DL.test(y,B,p)
}
\arguments{
  \item{y}{ financial return time series }
  \item{B}{ the number of bootstrap iterations, the default is 300}
  \item{p}{ the lag value, the default is 1}
}

\value{
  \item{Cp}{Cramer von Mises test statistic}
  \item{Kp}{Kolmogorov-Smirnov test statistic}
  \item{Cp_pval}{wild bootstrap p-value of the Cp test}
  \item{Kp_pval}{wild bootstrap p-value of the Kp test}
}
\references{ 
Domingues M.A. and Lobato, I. N., 2003,  Testing the Martingale Difference Hypothesis, Econometrics Reviews, 22, p351-377.

Charles, A. Darne, O. Kim, J.H. 2011, Small Sample Proeprties of Alternative Tests for Martingale Difference Hypothesis, Economics Letters, 110(2), 151-154.}
\author{ Jae H. Kim}

\examples{
r <- rnorm(50)           
DL.test(r,B=100)
# B=100 is used for fast execution in the example. 
# Use a higher number in actual application
}
\keyword{  htest  }
