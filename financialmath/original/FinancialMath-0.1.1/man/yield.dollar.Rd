\arguments{
  \item{cf}{vector of cash flows }
  \item{times}{vector of times for when cash flows occur}
  \item{start}{beginning balance}
  \item{end}{ending balance}
  \item{endtime}{end time of comparison}
}
\description{Calculates the dollar weighted yield.}
\name{yield.dollar}
\alias{yield.dollar}
\title{Dollar Weighted Yield}
\keyword{yield}
\note{Time of comparison (endtime) must be larger than any number in vector of cash flow times.

Length of cashflow vector and times vector must be equal.}
\value{The dollar weighted yield.}
\examples{yield.dollar(cf=c(20,10,50),times=c(.25,.5,.75),start=100,end=175,endtime=1)

yield.dollar(cf=c(500,-1000),times=c(3/12,18/12),start=25200,end=25900,endtime=21/12)}
\usage{yield.dollar(cf, times, start, end, endtime)}
\seealso{
  \code{\link{yield.time}}
}
\details{
\eqn{I=end-start-\sum_{k=1}^ncf_k}

\eqn{i^{dw}=\frac{I}{start*endtime-\sum_{k=1}^ncf_k*(endtime-times_k)}}
}
