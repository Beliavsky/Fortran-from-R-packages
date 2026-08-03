\name{maximType}
\alias{maximType}
\alias{maximType.default}
\alias{maximType.maxim}
\alias{maximType.MLEstimate}
\title{Type of Minimization/Maximization}
\description{
  Returns the type of optimization as supplied by the optimisation routine.
}
\usage{
  maximType(x)
}
\arguments{
  \item{x}{object of class 'maxim' or another object which
    involves numerical optimisation.
  }
}
\value{
  A text message, describing the involved optimisation algorithm
}
\author{Ott Toomet}
\seealso{\code{\link{maxNR}}}
\examples{
## maximize two-dimensional exponential hat.  True maximum c(2,1):
f <- function(a) exp(-(a[1] - 2)^2 - (a[2] - 1)^2)
m <- maxNR(f, start=c(0,0))
coef(m)
maximType(m)
## Now use BFGS maximisation.
m <- maxBFGS(f, start=c(0,0))
maximType(m)
}
\keyword{optimize}
\keyword{methods}
