\name{symm.huber.wt}
\alias{symm.huber.wt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Weighted Symmetrized Huber Scatter Matrix }
\description{
  Iterative algorithm to estimate the weighted symmetrized Huber scatter matrix.
}
\usage{
symm.huber.wt(X, wt = rep(1, nrow(X)), qg = 0.9, init = NULL, 
              eps = 1e-06, maxiter = 100, na.action = na.fail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ numeric data frame or matrix.}
  \item{wt}{vector of weights. Should be nonnegative and at least one larger than zero.}
  \item{qg}{tuning parameter. Should be between 0 and 1. The default is 0.9.}
  \item{init}{an optional matrix giving the starting value for the iteration.}
  \item{eps}{convergence tolerance.}
  \item{maxiter}{maximum number of iterations.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}
\details{
  The weighted symmetrized Huber scatter matrix is the regular Huber scatter matrix for the weighted pairwise differences of the observations taken wrt to the origin.

  Note that this function is memory comsuming and slow for large data sets since the matrix is based on all pairwise difference of the observations.
   
}
\value{
  a matrix.
}
\references{ \cite{Sirkia, S., Taskinen, S. and Oja, H. (2007), Symmetrised M-estimators of scatter. \emph{Journal of Multivariate Analysis}, \bold{98}, 1611--1629. } }
\author{Klaus Nordhausen}

\seealso{ \code{\link{symm.huber}} }
\examples{
set.seed(1)
cov.matrix.1 <- matrix(c(3,2,1,2,4,-0.5,1,-0.5,2), ncol = 3)
X.1 <- rmvnorm(100, c(0,0,0), cov.matrix.1)
cov.matrix.2 <- diag(1,3)
X.2 <- rmvnorm(50, c(1,1,1), cov.matrix.2)
X <- rbind(X.1, X.2)

D1 <-  symm.huber.wt(X, rep(c(0,1), c(100,50)))
D2 <-  symm.huber.wt(X, rep(c(1,0), c(100,50)))

D1
D2

rm(.Random.seed)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{robust}
