\name{compoundInterest}
\alias{compoundInterest}
\alias{simple2logReturns}
\title{Compute compound interest}
\description{
  Compute compound interest for a given number of
  periods, compounding with an indicated frequency per period.  
}
\usage{
compoundInterest(interest, periods = 1, frequency = 1, net.value = FALSE)

simple2logReturns(R)
}
\arguments{
  \item{interest}{ rate of interest per period (usually per year). }
  \item{periods}{ number of periods over which to compound. }
  \item{frequency}{ number of times per period to compound;
         \code{frequency = Inf} to convert simple to log returns.}
  \item{net.value}{ if TRUE, return the total value per
    unit invested;  otherwise return net increase =
    (net value - 1). }
  \item{R}{ simple interest to be converted to log(returns). }
}
\details{
  
  These functions are vectorized for all arguments.  The code uses
  optionally \code{expm1(x) = exp(x) - 1} and \code{log1p(x) = log(1+x)}
  which can preserve numerical precision for x very close to 0.

}
\value{
  vector of the length of the longest argument.  
}
\references{
  Ruey Tsay (2005)
  Analysis of Financial Time Series, 2nd ed. (Wiley, p. 6)
}
%\seealso{ \code{\link{}} }
\examples{
# "Net Value" column of Tsay Table 1.1, p. 4
compoundInterest(0.1, frequency = c(1, 2, 4, 12, 52, 365, Inf),
                      net.value = FALSE)
# Example 1.1, p. 6
compoundInterest(.0446, freq = Inf)
# Inverse of Example 1.1
simple2logReturns(.0456)
}
\keyword{ arith }
