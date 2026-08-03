\encoding{latin1}

\name{tyler.shape}
\alias{tyler.shape}

\title{ Tyler's Shape Matrix}
\description{
  Iterative algorithm to estimate Tyler's shape matrix.
}
\usage{
tyler.shape(X, location = NULL, init = NULL, steps = Inf, eps = 1e-06, 
            maxiter = 100, in.R = FALSE, print.it = FALSE, 
            na.action = na.fail)
}

\arguments{
  \item{X}{numeric data matrix or dataframe.}
  \item{location}{if NULL the sample mean is used, otherwise a vector with the location can be specified.}
  \item{init}{an optional matrix giving the starting value for the iteration}
  \item{steps}{a fixed number of iteration steps to take. See details.}
  \item{eps}{convergence tolerance.}
  \item{maxiter}{maximum number of iterations.}
  \item{in.R}{logical. If TRUE R-code (and not C) is used in the iteration}
  \item{print.it}{logical. If TRUE prints the number of iterations, otherwise not.}
  \item{na.action}{ a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}
\details{The most robust M-estimator of shape. It is proportional to the regular covariance matrix for elliptical contoured distributions. 
The estimate is in such a way standardized, that its determinate is 1.

The algorithm requires an estimate of location, if none is provided, the sample mean is used. Observations which are equal to the location estimate are removed form the data.

The algorithm also allows for a k-step version where the iteration is run for a fixed number of steps instead of until convergence. If \code{steps} is finite that number of steps is taken and \code{maxiter} is ignored.

A different implementation is available in the package fastM as the function \code{TYLERshape}.
}
\value{
  A matrix.
}
\references{
   \cite{Tyler, D.E. (1987), A distribution-free M-estimator of scatter, \emph{Annals of Statistics}, \bold{15},  234--251. }  
}
\author{Klaus Nordhausen, and Seija Sirkia}
\seealso{ \code{\link{duembgen.shape}}, \code{\link{HR.Mest}}}
\examples{
set.seed(654321)
cov.matrix <- matrix(c(3,2,1,2,4,-0.5,1,-0.5,2), ncol=3)
X <- rmvnorm(100, c(0,0,0), cov.matrix)
tyler.shape(X)
tyler.shape(X, location=0)
cov.matrix/det(cov.matrix)^(1/3)
rm(.Random.seed)
}
\keyword{multivariate}
\keyword{robust}

