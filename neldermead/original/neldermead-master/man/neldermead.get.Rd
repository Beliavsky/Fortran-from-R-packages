\name{neldermead.get}
\alias{neldermead.get}
\title{
  Get the value for the given element
}
\description{
  Get the value for the given element in a neldermead object.
}
\usage{
  neldermead.get(this = NULL, key = NULL)
}
\arguments{
  \item{this}{A neldermead object.}
  \item{key}{The name of the key to query.}
}
\value{
  Return the value of the list element \code{key}, or an error message if
  \code{key} does not exist in the neldermead object \code{this}.
}
\author{
Author of Scilab neldermead module: Michael Baudin (INRIA - Digiteo)

Author of R adaptation: Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}
\seealso{
  \code{\link[neldermead]{neldermead.set}},
  \code{\link[optimbase]{optimbase.get}}
}
\keyword{ method }
