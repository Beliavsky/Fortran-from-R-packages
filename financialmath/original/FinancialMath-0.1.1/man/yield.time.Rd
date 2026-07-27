\arguments{
  \item{cf}{vector of cash flows}
  \item{bal}{vector of balances}
}
\name{yield.time}
\alias{yield.time}
\title{Time Weighted Yield}
\description{Calculates the time weighted yield.}
\usage{yield.time(cf,bal)}
\author{Kameron Penn and Jack Schmidt}
\keyword{yield}
\keyword{time}
\note{Length of cash flows must be one less than the length of balances.

If lengths are equal, it will not use final cash flow.}
\value{The time weighted yield.}
\examples{yield.time(cf=c(0,200,100,50),bal=c(1000,800,1150,1550,1700))}
\seealso{
  \code{\link{yield.dollar}}
}
\details{
\eqn{i^{tw}=\prod_{k=1}^n (\frac{bal_{1+k}}{bal_k+cf_k})-1}
}
