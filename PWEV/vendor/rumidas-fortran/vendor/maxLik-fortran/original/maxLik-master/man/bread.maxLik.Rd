\name{bread.maxLik}
\alias{bread}
\alias{bread.maxLik}

\title{Bread for Sandwich Estimator}
\description{
   Extracting an estimator for the \sQuote{bread} of the sandwich estimator,
   see \code{\link[sandwich]{bread}}.
}
\usage{
\method{bread}{maxLik}( x, ... )
}
\arguments{
   \item{x}{an object of class \code{maxLik}.}
   \item{\dots}{further arguments (currently ignored).}
}
\value{
   Matrix, the inverse of the expectation of the second derivative
   (Hessian matrix) of the log-likelihood function
   with respect to the parameters.
   In case of the simple Maximum Likelihood, it is
   equal to the variance covariance matrix of the parameters,
   multiplied by the number of observations.
}
\section{Warnings}{
  The \pkg{sandwich} package is required for this function. 

  This method works only if the observaton-specific gradient information
  was available for the estimation.  This is the case if the
  observation-specific gradient was supplied (see the \code{grad}
  argument for \code{\link{maxLik}}), or the log-likelihood function
  returns a vector of observation-specific values.
}
\author{
   Arne Henningsen
}

\seealso{\code{\link[sandwich]{bread}}, \code{\link{maxLik}}.}

\examples{
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t

## Estimate with numeric gradient and hessian
a <- maxLik(loglik, start=1 )

# Extract the "bread"
library( sandwich )
bread( a )

all.equal( bread( a ), vcov( a ) * nObs( a ) )
}

\keyword{methods}
