\name{summary.onls}
\alias{summary.onls}

\title{Summary function for 'onls' objects}

\description{
Provides a summary for the parameters of the converged fit including their error estimates, and the standard errors in respect to vertical and orthogonal residuals.
}

\usage{
\method{summary}{onls}(object, correlation = FALSE, symbolic.cor = FALSE, ...)
}

\arguments{
  \item{object}{an object returned from \code{\link{onls}}.}
  \item{correlation}{logical. If \code{TRUE}, the correlation matrix of the estimated parameters is returned and printed.}
  \item{symbolic.cor}{logical. If \code{TRUE}, print the correlations in a symbolic form.}
  \item{...}{further arguments passed to or from other methods.}       
}

\value{
A summary similar to \code{\link{summary.nls}} containing the formula, parameter estimates with their errors/significance, and standard errors in respect to vertical and orthogonal residuals.
}
     
\author{
Andrej-Nikolai Spiess
}
        
\examples{
## See 'onls'.
}

\keyword{models}
\keyword{nonlinear}
