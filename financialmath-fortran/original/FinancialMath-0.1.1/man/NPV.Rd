\arguments{
  \item{cf0}{cash flow at period 0}
  \item{cf}{vector of cash flows}
  \item{times}{vector of the times for each cash flow}
  \item{i}{interest rate per period}
  \item{plot}{tells whether or not to plot the time diagram of the cash flows}
}
\description{Calculates the net present value for a series of cash flows, and provides a time diagram of the cash flows.}
\name{NPV}
\alias{NPV}
\title{Net Present Value}
\keyword{value}
\note{The periods in t must be positive integers.

The lengths of cf and t must be equal.}
\value{The NPV.}
\examples{NPV(cf0=100,cf=c(50,40),times=c(3,5),i=.01)

NPV(cf0=100,cf=50,times=3,i=.05)

NPV(cf0=100,cf=c(50,60,10,20),times=c(1,5,9,9),i=.045)}
\usage{NPV(cf0,cf,times,i,plot=FALSE)}
\seealso{
  \code{\link{IRR}}
}
\details{
\eqn{NPV=cf0-\sum_{k=1}^n\frac{cf_k}{(1+i)^{times_k}}}
}
