\name{fminbnd.outputfun}
\alias{fminbnd.outputfun}
\title{
  fminbnd Output Function Call
}
\description{
  This function calls the output function and make it match neldermead
  requirements. It is used in the \code{fminbnd} function as the
  \code{outputcommand} element of the neldermead object (see
  \code{?neldermead} and \code{?neldermead.set}).
}
\usage{
  fminbnd.outputfun(state = NULL, data = NULL, fmsdata = NULL)
}
\arguments{
  \item{state}{The current state of the algorithm either 'init', 'iter' or
    'done'.}
  \item{data}{The data at the current state. This is an object of class 
    'neldermead.data', i.e. a list with the following elements: \describe{
      \item{x}{The current parameter estimates.}
      \item{fval}{The current value of the cost function.}
      \item{simplex}{The current simplex object.}
      \item{iteration}{The number of iterations performed.}
      \item{funccount}{The number of function evaluations.}
      \item{step}{The type of step in the previous iteration.}
    }
  }
  \item{fmsdata}{This is an object of class 'optimbase.functionargs' which 
    contains specific data of the \code{fminbnd} algorithm: \describe{
      \item{Display}{what to display}
      \item{OutputFcn}{the array of output functions}
      \item{PlotFcns}{the array of plot functions}
    }
  }
}
\value{
  This function does not return any data, but execute the output function(s).
}
\author{
Author of Scilab neldermead module: Michael Baudin (INRIA - Digiteo)

Author of R adaptation: Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}
\seealso{
  \code{\link[neldermead]{fminbnd}},
  \code{\link[neldermead]{neldermead}},
  \code{\link[neldermead]{neldermead.set}},
}
\keyword{ method }
