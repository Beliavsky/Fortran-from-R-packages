\name{duembgen.shape}
\alias{duembgen.shape}

\title{ Duembgen's Shape Matrix }
\description{
  Iterative algorithm to estimate Duembgen's shape matrix. 
}
\usage{
duembgen.shape(X, init = NULL, steps = Inf, eps = 1e-06, 
               maxiter = 100, in.R = FALSE, na.action = na.fail, ...)
}

\arguments{
  \item{X}{numeric data matrix or dataframe.}
  \item{init}{an optional matrix giving the starting value for the iteration. Otherwise the regular covariance is used after transforming it to a shape matrix wit determinant 1.}
  \item{steps}{a fixed number of iteration steps to take. See details.}
  \item{eps}{ convergence tolerance.}
  \item{maxiter}{ maximum number of iterations.}
  \item{in.R}{logical. If TRUE R-code (and not C) is used in the iteration}
  \item{na.action}{ a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
  \item{\dots}{ other arguments passed on to  \code{\link{tyler.shape}}.}
}
\details{Duembgen's shape matrix can be seen as \code{\link{tyler.shape}}'s matrix wrt to the origin for the pairwise differences of the observations.
Therefore this shape matrix needs no location parameter.

The function is, however, slow if the dataset is large.

The algorithm also allows for a k-step version where the iteration is run for a fixed number of steps instead of until convergence. If \code{steps} is finite that number of steps is taken and \code{maxiter} is ignored.

A better implementation is available in the package fastM as the function \code{DUEMBGENshape}.
}
\value{
  A matrix.
}
\references{
   \cite{Duembgen, L. (1998), On Tyler's M-functional of scatter in high dimension, \emph{Annals of Institute of Statistical Mathematics}, \bold{50},  471--491. }  
}
\author{Klaus Nordhausen, Seija Sirkia, and some of the C++ is based on work by Jari Miettinen}
\seealso{ \code{\link{tyler.shape}}, \code{\link{duembgen.shape.wt}}}
\examples{
set.seed(654321)
cov.matrix <- matrix(c(3,2,1,2,4,-0.5,1,-0.5,2), ncol=3)
X <- rmvnorm(100, c(0,0,0), cov.matrix)
cov.matrix/det(cov.matrix)^(1/3)
duembgen.shape(X)
rm(.Random.seed)
}
\keyword{multivariate}
\keyword{robust}
