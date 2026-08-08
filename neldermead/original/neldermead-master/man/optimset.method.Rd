\name{optimset.method}
\alias{optimset.method}
\title{
  Default set of optimization options
}
\description{
  This function returns a default set of optimization options for defined
  'methods'; \code{optimset.method} is called by \code{optimset} when a
  \code{method} was provided as input. Currently, the only valid \code{method}
  is 'fminsearch'.
}
\usage{
  optimset.method(method = NULL)
}
\arguments{
  \item{method}{A character string.}
}
\value{
  Returns a list with the following fields: Display, FunValCheck, MaxFunEvals,
  MaxIter, OutputFcn, PlotFcns, TolFun, and TolX.
}
\author{
  Author of Scilab neldermead module: Michael Baudin (INRIA - Digiteo)

  Author of R adaptation: Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}
\seealso{
  \code{\link{optimset}}
}
\examples{
  optimset.method('fminsearch')
  # Will fail
  try(optimset.method('abc'))
}
\keyword{ method }
