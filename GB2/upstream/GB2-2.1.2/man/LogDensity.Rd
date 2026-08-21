\name{LogDensity}
\Rdversion{2.1}
\alias{LogDensity}
\alias{logf.gb2}
\alias{dlogf.gb2}
\alias{d2logf.gb2}
\title{
Log Density of the GB2 Distribution
}
\description{
Calculates the log density of the GB2 distribution for a single value or a vector of values. 
Calculates the first- and second-order partial derivatives of the log density evaluated at a single value.}
\usage{
logf.gb2(x, shape1, scale, shape2, shape3)
dlogf.gb2(xi, shape1, scale, shape2, shape3)
d2logf.gb2(xi, shape1, scale, shape2, shape3)
}
\arguments{
  \item{xi}{numeric; a data value.}
  \item{x}{numeric; a vector of data values.}
  \item{shape1}{numeric; positive parameter.}
  \item{scale}{numeric; positive parameter.}
  \item{shape2, shape3}{numeric; positive parameters of the Beta distribution.}
}
\details{
We calculate \eqn{log(f(x, \theta))}, where \eqn{f} is the GB2 density with parameters \code{shape1} \eqn{= a}, \code{scale} \eqn{= b}, 
\code{shape2} \eqn{= p} and \code{shape3} \eqn{= q}, \eqn{\theta} is the parameter vector. We calculate the first- and second-order partial derivatives of \eqn{log(f(x, \theta))} with 
respect to the parameter vector \eqn{\theta}.
}
\value{
Depending on the input \code{logf.gb2} gives the log density for a single value or a vector of values. \code{dlogf.gb2} gives the vector of the four first-order partial derivatives of the log density and 
\code{d2logf.gb2} gives the \eqn{4 \times 4} matrix of second-order partial derivatives of the log density.
}
\references{
Brazauskas, V. (2002)
 Fisher information matrix for the Feller-Pareto distribution.
 \emph{Statistics & Probability Letters}, \bold{59}, 159--167.
}
\author{
Desislava Nedyalkova
}
\keyword{distribution}

