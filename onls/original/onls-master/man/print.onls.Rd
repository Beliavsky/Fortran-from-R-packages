\name{print.onls}
\alias{print.onls}

\title{Printing function for 'onls' objects}

\description{
Provides a printed summary of the converged fit obtained from \code{\link{onls}}.
}

\usage{
\method{print}{onls}(x, ...)
}

\arguments{
  \item{x}{an object returned from \code{\link{onls}}.}
  \item{...}{other parameters for future methods.}       
}

\value{
A printed summary containing the formula, data name, converged parameters, vertical residual sum-of-squares \deqn{\sum_{i=1}^{n}(y_i - \hat{y_i})^2}, orthogonal residual sum-of-squares \deqn{\sum_{i=1}^{n}\| D_i \|^2}, number of points with orthogonality as obtained from \code{\link{check_o}} and number of iterations to convergence. 
}
     
\author{
Andrej-Nikolai Spiess
}
        
\examples{
## See 'onls'.
}

\keyword{models}
\keyword{nonlinear}
