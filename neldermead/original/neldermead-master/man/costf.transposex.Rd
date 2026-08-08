\name{costf.transposex}
\alias{costf.transposex}
\title{
  Cost Function Call
}
\description{
  Call the cost function after transposition of the value of the point
  estimate x, so that the input row vector, given by optimsimplex, is
  transposed into a  column vector as required by the cost function.
}
\usage{
  costf.transposex(x = NULL, this = NULL)
}
\arguments{
  \item{x}{The point estimate provide as a row matrix.}
  \item{this}{A neldermead object.}
}
\value{
  Return the value of the cost function (called by \code{neldermead.costf}.)
}
\author{
Author of Scilab neldermead module: Michael Baudin (INRIA - Digiteo)

Author of R adaptation: Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}
\seealso{
  \code{\link[neldermead]{neldermead.costf}}
}
\keyword{ method }
