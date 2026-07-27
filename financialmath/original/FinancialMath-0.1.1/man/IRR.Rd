\arguments{
  \item{cf0}{cash flow at period 0}
  \item{cf}{vector of cash flows}
  \item{times}{vector of the times for each cash flow}
  \item{plot}{option whether or not to provide the time diagram}
}
\name{IRR}
\alias{IRR}
\title{Internal Rate of Return}
\description{Calculates internal rate of return for a series of cash flows, and provides a time diagram of the cash flows.}
\usage{IRR(cf0,cf,times,plot=FALSE)}
\author{Kameron Penn and Jack Schmidt}
\keyword{irr}
\keyword{yield}
\note{Periods in t must be positive integers.

Uses polyroot function to solve equation given by series of cash flows, meaning that in the case of having a negative IRR, multiple answers may be returned.
}
\value{The internal rate of return.}
\examples{IRR(cf0=1,cf=c(1,2,1),times=c(1,3,4))

IRR(cf0=100,cf=c(1,1,30,40,50,1),times=c(1,1,3,4,5,6))}
\seealso{
  \code{\link{NPV}}
}
\details{
\eqn{cf0=\sum_{k=1}^n\frac{cf_k}{(1+irr)^{times_k}}}
}
