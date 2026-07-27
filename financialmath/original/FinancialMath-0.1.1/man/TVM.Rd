\arguments{
  \item{pv}{present value}
  \item{fv}{future value}
  \item{n}{number of periods}
  \item{i}{nominal interest rate convertible ic times per period}
  \item{ic}{interest conversion frequency per period}
  \item{plot}{tells whether or not to produce a plot of the time value at each period}
}
\description{Solves for the present value, future value, time, or the interest rate for the accumulation of money earning compound interest.  It can also plot the time value for each period.}
\name{TVM}
\alias{TVM}
\title{Time Value of Money}
\keyword{time}
\keyword{value}
\note{Exactly one of pv, fv, n, or i must be NA (unknown).}
\value{Returns a matrix of the input variables and calculated unknown variables.}
\examples{
TVM(pv=10,fv=20,i=.05,ic=2,plot=TRUE)

TVM(pv=50,n=5,i=.04,plot=TRUE)
}
\usage{TVM(pv=NA,fv=NA,n=NA,i=NA,ic=1,plot=FALSE)}
\seealso{
  \code{\link{cf.analysis}}
}
\details{
\eqn{j=(1+\frac{i}{ic})^{ic}-1}

\eqn{fv=pv*(1+j)^n}
}
