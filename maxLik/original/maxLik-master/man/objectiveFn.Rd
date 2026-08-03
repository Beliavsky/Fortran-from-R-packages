\name{objectiveFn}
\alias{objectiveFn}
\alias{objectiveFn.maxim}
\title{Optimization Objective Function}
\description{
  This function returns the optimization objective function from a
  \sQuote{maxim} object.
}
\usage{
objectiveFn(x, \dots)
\method{objectiveFn}{maxim}(x, \dots)
}
\arguments{
  \item{x}{an optimization result, inheriting from class \sQuote{maxim}}
  \item{\dots}{other arguments for methods}
}
\value{
  function, the function that was optimized.  It can be directly called,
  given that all necessary variables are accessible from the current environment.
}
\author{Ott Toomet}
\examples{
hatf <- function(theta) exp(- theta \%*\% theta)
res <- maxNR(hatf, start=c(0,0))
print(summary(res))
print(objectiveFn(res))
print(objectiveFn(res)(2)) # 0.01832
}
\keyword{methods}
\keyword{optimize}
