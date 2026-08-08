\name{neldermead.function}
\alias{neldermead.function}
\title{
  Call Cost Function.
}
\description{
  Simple way to compute the value of the cost function specified in a neldermead
  object.
}
\usage{
  neldermead.function(this = NULL, x = NULL)
}
\arguments{
  \item{this}{A neldermead object.}
  \item{x}{The point estimate where the cost function is to be evaluated.}
}
\value{
  Returns the value of the cost function.
}
\author{
Author of Scilab neldermead module: Michael Baudin (INRIA - Digiteo)

Author of R adaptation: Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}
\keyword{ method }
