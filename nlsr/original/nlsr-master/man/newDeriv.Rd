\name{newDeriv}
\alias{newDeriv}
\title{newDeriv}
\usage{
newDeriv (expr, deriv, derivEnv = sysDerivs)
}
\arguments{
\item{expr}{An expression represented as a quoted expression.}

\item{deriv}{The derivative.  See Details below. }

\item{derivEnv}{The environment in which to store the setting.}
}
\description{
Define a new derivative expression.
}
\details{
Both \code{expr} and \code{deriv} are treated
as unevaluated expressions.

\code{deriv} should include the sum of derivatives with respect to
all variables similar to a \dQuote{total derivative} using
\code{\link{D}} for the differential.  See the examples
below.

This is mainly intended for internal use.
}
\examples{
newDeriv(sin(x), cos(x)*D(x))
newDeriv(x*y, x*D(y) + D(x)*y)

nlsDeriv(quote(sin(x)*sin(y)), "x")
}
\keyword{internal}