\name{returnCode}
\alias{returnCode}
\alias{returnCode.default}
\alias{returnCode.maxLik}
\alias{returnMessage}
\alias{returnMessage.default}
\alias{returnMessage.maxim}
\alias{returnMessage.maxLik}
\title{Success or failure of the optimization}
\description{
  These function extract success or failure information from optimization objects.
  The \code{returnCode} gives a numeric code, and \code{returnMessage} a
  brief description about the success or
  failure of the optimization, and point to the problems occured (see
  documentation for the
  corresponding functions).
}
\usage{
returnCode(x, ...)
\method{returnCode}{default}(x, ...)
\method{returnCode}{maxLik}(x, ...)
returnMessage(x, ...)
\method{returnMessage}{maxim}(x, ...)
\method{returnMessage}{maxLik}(x, ...)
}
\arguments{
  \item{x}{object, usually an optimization result}
  \item{...}{further arguments for other methods}
}
\details{
  \code{returnMessage} and \code{returnCode} are a generic functions, with methods
  for various optimisation algorithms.
  The message should either describe
  the convergence (stopping condition),
  or the problem.

  The known codes and the related messages are:
  \describe{
    \item{1}{ gradient close to zero (normal convergence).}
    \item{2}{ successive function values within tolerance limit (normal
      convergence).}
    \item{3}{ last step could not find higher value (probably not
      converged).  This is related to line search step getting too
      small, usually because hitting the boundary of the parameter
      space.  It may also be related to attempts to move to a wrong
      direction because of numerical errors.  In some cases it can be
      helped by changing \code{steptol}.}
    \item{4}{ iteration limit exceeded.}
    \item{5}{ Infinite value.}
    \item{6}{ Infinite gradient.}
    \item{7}{ Infinite Hessian.}
    \item{8}{Successive function values withing relative tolerance
      limit (normal convergence).}
    \item{9}{ (BFGS) Hessian approximation cannot be improved because of
      gradient did not change.  May be related to numerical
      approximation problems or wrong analytic gradient.
    }
    \item{10}{
      Lost patience: the optimizer has hit an inferior value too many
      times (see \code{\link{maxSGA}} for more information)
    }
    \item{100}{ Initial value out of range.}
    }
}
\value{
  Integer for \code{returnCode}, character for \code{returnMessage}.
  Different optimization routines may define it in a different way.
}
\author{Ott Toomet}
\seealso{\code{\link{maxNR}}, \code{\link{maxBFGS}}}
\examples{
## maximise the exponential bell
f1 <- function(x) exp(-x^2)
a <- maxNR(f1, start=2)
returnCode(a) # should be success (1 or 2)
returnMessage(a)
## Now try to maximise log() function
a <- maxNR(log, start=2)
returnCode(a) # should give a failure (4)
returnMessage(a)
}
\keyword{methods}
\keyword{utilities}
