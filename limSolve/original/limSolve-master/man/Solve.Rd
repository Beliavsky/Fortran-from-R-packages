\name{Solve}
\alias{Solve}
\title{
  Generalised inverse solution of Ax = B
}
\description{
  Generalised inverse solution of \deqn{Ax=B}
}

\note{  
\code{Solve} uses the Moore-Penrose generalized inverse of matrix \code{A}
  (function \code{\link[MASS]{ginv}} from package MASS).
 
  \code{\link{solve}}, the \R default requires a square, positive
  definite A. \code{Solve} does not have this restriction.
}
\usage{
Solve (A, B = diag(nrow = nrow(A)), tol = sqrt(.Machine$double.eps))
}

\arguments{
  \item{A }{numeric matrix containing the coefficients of the equations
    \eqn{Ax=B}.
  }
  \item{B }{numeric matrix containing the right-hand sides of the equations;
    the default is the unity matrix, in which case the function will return
    the Moore-Penrose generalized inverse of matrix A.
  }
  \item{tol }{tolerance for selecting singular values.
  }
}

\value{
  a vector with the generalised inverse solution.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
A <- matrix(nrow = 4, ncol = 3, data = c(1:8, 6, 8, 10, 12)) # col3 = col1+col2
B <- 0:3
X <- Solve(A, B)            # generalised inverse solution
A \%*\% X - B                 # should be zero (except for roundoff)
(gA <- Solve(A))           # generalised inverse of A
}

\seealso{
  \code{\link[MASS]{ginv}} to estimate the Moore-Penrose generalized inverse
  of a matrix, in package \code{MASS},

  \code{\link{solve}} the R default
}
\references{
package \code{MASS}:

  Venables, W. N. & Ripley, B. D. (2002) Modern Applied
  Statistics with S. Fourth Edition. Springer, New
  York. ISBN 0-387-95457-0
}
\keyword{array}

