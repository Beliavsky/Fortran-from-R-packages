\name{gradient}
\alias{gradient}
\alias{gradient.maxim}
\alias{estfun}
\alias{estfun.maxLik}

\title{Extract Gradients Evaluated at each Observation}
\description{
   Extract the gradients of the log-likelihood function evaluated
   at each observation (\sQuote{Empirical Estimating Function},
   see \code{\link[sandwich]{estfun}}).
}
\usage{
\method{estfun}{maxLik}(x, ...)
\method{gradient}{maxim}(x, ...)
}
\arguments{
   \item{x}{an object inheriting from class \code{maxim} (for \code{gradient})
     or \code{maxLik}. (for \code{estfun}.)}
   \item{\dots}{further arguments (currently ignored).}
}
\value{
  \item{\code{gradient}}{vector, objective function gradient at estimated
    maximum (or the last calculated value
    if the estimation did not converge.)}
  \item{\code{estfun}}{
    matrix, observation-wise log-likelihood gradients
   at the estimated parameter value
   evaluated at each observation.  Observations in rows, parameters in columns.}
}
\section{Warnings}{
   The \pkg{sandwich} package must be loaded in order to use \code{estfun}.

   \code{estfun} only works if the observaton-specific gradient information
   was available for the estimation.  This is the case of the
   observation-specific gradient was supplied (see the \code{grad}
   argument for \code{\link{maxLik}}), or the log-likelihood function
   returns a vector of observation-specific values.
}
\author{
   Arne Henningsen, Ott Toomet
}

\seealso{\code{\link{hessian}}, \code{\link[sandwich]{estfun}}, \code{\link{maxLik}}.}

\examples{
## ML estimation of exponential duration model:
t <- rexp(10, 2)
loglik <- function(theta) log(theta) - theta*t

## Estimate with numeric gradient and hessian
a <- maxLik(loglik, start=1 )

gradient(a)
# Extract the gradients evaluated at each observation
library( sandwich )
estfun( a )

## Estimate with analytic gradient.
## Note: it returns a vector
gradlik <- function(theta) 1/theta - t
b <- maxLik(loglik, gradlik, start=1)
gradient(a)
estfun( b )
}

\keyword{methods}
