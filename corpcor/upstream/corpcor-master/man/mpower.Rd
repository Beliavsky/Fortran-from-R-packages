\name{mpower}
\alias{mpower}

\title{Compute the Power of a Real Symmetric Matrix}
\usage{
mpower(m, alpha, pseudo=FALSE, tol)
}
\arguments{
  \item{m}{a real-valued symmetric matrix.}
  
  \item{alpha}{exponent.}

  \item{pseudo}{if \code{pseudo=TRUE} then all zero eigenvalues are dropped 
                (e.g. for computing the pseudoinverse).  The default is
                to use all eigenvalues.}

\item{tol}{tolerance - eigenvalues with absolute value smaller or equal
           to \code{tol} are considered identically zero (default:
	     \code{tol = max(dim(m))*max(abs(eval))*.Machine$double.eps}).
	     }

}
\description{
  \code{mpower} computes \eqn{m^alpha}, i.e.
  the \code{alpha}-th power of the real symmetric
   matrix \code{m}.
}

\value{
  \code{mpower} returns
  a matrix of the same dimensions as \code{m}.
}

\details{
   The matrix power of \code{m} is obtained by first computing the spectral
   decomposition of \code{m}, and subsequent modification of the resulting eigenvalues.

  Note that \code{m} is assumed to by symmetric, and only 
  its lower triangle (diagonal  included) is used in \code{\link{eigen}}.

  For computing the matrix power of \code{\link{cor.shrink}}  use
  the vastly more efficient function  \code{\link{powcor.shrink}}.
}

\author{
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
}
\seealso{\code{\link{powcor.shrink}}, \code{\link{eigen}}.}

\examples{
# load corpcor library
library("corpcor")

# generate symmetric matrix
p = 10
n = 20
X = matrix(rnorm(n*p), nrow = n, ncol = p)
m = cor(X)

m \%*\% m
mpower(m, 2)

solve(m)
mpower(m, -1)

msq = mpower(m, 0.5)
msq \%*\% msq
m

mpower(m, 1.234)
}
\keyword{algebra}
