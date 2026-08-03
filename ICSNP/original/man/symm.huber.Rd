\name{symm.huber}
\alias{symm.huber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Symmetrized Huber Scatter Matrix}
\description{
  Iterative algorithm to estimate the symmetrized Huber scatter matrix.
}
\usage{
symm.huber(X, qg = 0.9, init = NULL, eps = 1e-06, maxiter = 100, 
           na.action = na.fail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ numeric data frame or matrix.}
  \item{qg}{tuning parameter. Should be between 0 and 1. The default is 0.9.}
  \item{init}{an optional matrix giving the starting value for the iteration.}
  \item{eps}{convergence tolerance.}
  \item{maxiter}{maximum number of iterations.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}
\details{
  The symmetrized Huber scatter matrix is the regular Huber scatter matrix for the pairwise differences of the observations taken wrt to the origin.
  
  Note that this function might be memory comsuming and slow for large data sets since the matrix is based on all pairwise difference of the observations.
  
  The function \code{symmhuber} in the package \pkg{SpatialNP} offers also a k-step option. The \pkg{SpatialNP} package contains also the function \code{mvhuberM} for the regular multivariate Huber location
  and scatter estimatior.
}
\value{
  a matrix.
}
\references{ \cite{Sirkia, S., Taskinen, S. and Oja, H. (2007), Symmetrised M-estimators of scatter. \emph{Journal of Multivariate Analysis}, \bold{98}, 1611--1629. } }
\author{Klaus Nordhausen and Jari Miettinen}

\seealso{ \code{\link{symm.huber.wt}}, \code{symmhuber}, \code{mvhuberM}}
\examples{
set.seed(654321)
cov.matrix <- matrix(c(3,2,1,2,4,-0.5,1,-0.5,2), ncol=3)
X <- rmvnorm(100, c(0,0,0), cov.matrix)
symm.huber(X)
rm(.Random.seed)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{robust}
