\name{ddnorm}
\alias{ddnorm}

\title{Derivative of the Normal Distribution's Density Function}

\description{
   This function returns the derivative(s) of the density function
   of the normal (Gaussian) distribution with respect to the quantile,
   evaluated at the quantile(s), mean(s), and standard deviation(s)
   specified by arguments \code{x}, \code{mean}, and \code{sd}, respectively.
}

\usage{
ddnorm( x, mean = 0, sd = 1 )
}

\arguments{
  \item{x}{quantile or vector of quantiles.}
  \item{mean}{mean or vector of means.}
  \item{sd}{standard deviation or vector of standard deviations.}
}

\value{
   numeric value(s): derivative(s) of the density function
   of the normal distribution with respect to the quantile
}

\author{Arne Henningsen}

\seealso{\code{\link[stats]{dnorm}}}

\examples{
ddnorm( c( -1, 0, 1 ) )
}
\keyword{methods}
