\name{fminsearch.function}
\alias{fminsearch.function}
\title{
  fminsearch Cost Function Call
}
\description{
  This function calls the cost function and makes it match neldermead
  requirements. It is used in the \code{fminsearch} function as the
  \code{function} element of the neldermead object (see \code{?neldermead}
  and \code{?neldermead.set}). 
}
\usage{
  fminsearch.function(x = NULL, index = NULL, fmsfundata = NULL)
}
\arguments{
  \item{x}{A single column vector of parameter estimates.}
  \item{index}{An integer variable set to 2, indicating that only the cost
    function is to be computed by the algorithm.}
  \item{fmsfundata}{An object of class 'optimbase.functionargs' and with 
    (at least) a \code{fun} element, which contains the user-defined cost 
    function.}
}
\value{
  Returns a list with the following elements: \describe{
    \item{f}{The value of the cost function at the current point estimate.}
    \item{index}{The same \code{index} variable.}
    \item{this}{A list with a single element \code{costargument} which
      contains \code{fmsfundata}.}
  }
}
\author{
  Author of Scilab neldermead module: Michael Baudin (INRIA - Digiteo)

  Author of R adaptation: Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}
\seealso{
  \code{\link[neldermead]{fminsearch}},
  \code{\link[neldermead]{neldermead}},
  \code{\link[neldermead]{neldermead.set}},
}
\keyword{ method }
