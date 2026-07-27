\arguments{
  \item{cf}{vector of cash flows}
  \item{times}{vector of the periods for each cash flow}
  \item{i}{interest rate per period}
  \item{plot}{tells whether or not to plot the convexity}
  \item{time.d}{tells whether or not to plot the time diagram of the cash flows}
}
\description{Calculates the present value, macaulay duration and convexity, and modified duration and convexity for given cash flows.  It also plots the convexity and time diagram of the cash flows.}
\name{cf.analysis}
\alias{cf.analysis}
\title{Cash Flow Analysis}
\keyword{analysis}
\keyword{value}
\note{The periods in t must be positive integers.}
\value{A matrix of all of the calculated values.}
\examples{
cf.analysis(cf=c(1,1,101),times=c(1,2,3),i=.04,time.d=TRUE)

cf.analysis(cf=c(5,1,5,45,5),times=c(5,4,6,7,5),i=.06,plot=TRUE)}
\usage{cf.analysis(cf,times,i,plot=FALSE,time.d=FALSE)}
\seealso{
  \code{\link{TVM}}
}
\details{
\eqn{pv=\sum_{k=1}^n\frac{cf_k}{(1+i)^{times_k}}}

\eqn{MAC D=\frac{\sum_{k=1}^n times_k*(1+i)^{-times_k}*cf_k}{pv}}

\eqn{MOD D=\frac{\sum_{k=1}^n times_k*(1+i)^{-(times_k+1)}*cf_k}{pv}}

\eqn{MAC C=\frac{\sum_{k=1}^n {times_k}^2*(1+i)^{-times_k}*cf_k}{pv}}

\eqn{MOD C=\frac{\sum_{k=1}^n times_k*(times_k+1)*(1+i)^{-(times_k+2)}*cf_k}{pv}}

}
