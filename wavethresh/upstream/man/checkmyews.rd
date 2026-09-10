\name{checkmyews}
\alias{checkmyews}
\usage{
checkmyews(spec, nsim=10)
}
\arguments{
\item{spec}{The LSW spectrum}
\item{nsim}{The number of realizations}
}
\title{Check a LSW spectrum through repeated simulation and empirical averages}
\description{
Given a LSW spectrum this function simulates \code{nsim} realizations,
estimates the spectrum, and then averages the results. The large
sample averages should converge to the original spectrum.
}
\value{A LSW spectrum obtained as the average of \code{nsim} simulations
from the \code{spec} spectrum.}
\seealso{\code{\link{cns}},\code{\link{LSWsim}}, \code{\link{ewspec}}}
\keyword{ts}
\author{G P Nason}
