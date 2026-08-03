\name{vcov.maxLik}
\alias{vcov.maxLik}

\title{Variance Covariance Matrix of maxLik objects}
\description{
  Extract variance-covariance matrices from \code{\link{maxLik}} objects.
}
\usage{
   \method{vcov}{maxLik}( object, eigentol=1e-12, ... )
}
\arguments{
  \item{object}{a \sQuote{maxLik} object.}
  \item{eigentol}{
    eigenvalue tolerance, controlling when the Hessian matrix is
    treated as numerically singular.
  }
  \item{\dots}{further arguments (currently ignored).}
}
\value{
  the estimated variance covariance matrix of the coefficients.  In
  case of the estimated Hessian is singular, it's values are
  \code{Inf}.  The values corresponding to fixed parameters are zero.
}
\details{
  The standard errors are only calculated if the ratio of the smallest
  and largest eigenvalue of the Hessian matrix is less than
  \dQuote{eigentol}.  Otherwise the Hessian is treated as singular.
}
\author{
   Arne Henningsen,
   Ott Toomet
}

\seealso{\code{\link[stats]{vcov}}, \code{\link{maxLik}}.}

\examples{
## ML estimation of exponential random variables
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t
gradlik <- function(theta) 1/theta - t
hesslik <- function(theta) -100/theta^2
## Estimate with numeric gradient and hessian
a <- maxLik(loglik, start=1, control=list(printLevel=2))
vcov(a)
## Estimate with analytic gradient and hessian
a <- maxLik(loglik, gradlik, hesslik, start=1)
vcov(a)
}

\keyword{methods}
