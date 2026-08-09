\name{Solve.tridiag}
\alias{Solve.tridiag}
\title{
  Solution of a tridiagonal system of linear equations
}

\description{
  Solves the linear system of equations \deqn{Ax=B}
  
  where A has to be square and \emph{tridiagonal}, i.e with nonzero elements
  only on, one band above, and one band below the diagonal.
}

\usage{
Solve.tridiag ( diam1, dia, diap1, B=rep(0,times=length(dia)))
}

\arguments{
  \item{diam1 }{a vector with (nonzero) elements below the diagonal.
  }
  \item{dia }{a vector with (nonzero) elements on the diagonal.
  }
  \item{diap1 }{a vector with (nonzero) elements above the diagonal.
  }
  \item{B }{Right-hand side of the equations, a vector with length = number
   of rows of A, or a matrix with number of rows = number of rows of A.
  }
}
\value{
  matrix with the solution, \code{X}, of the tridiagonal system of equations Ax=B.
  The number of columns of this matrix equals the number of columns of B.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{# create tridagonal system: bands on diagonal, above and below
nn   <- 20                          # nr rows and columns of A
aa   <- runif(nn)
bb   <- runif(nn)
cc   <- runif(nn)

# full matrix
A                        <- matrix(nrow = nn, ncol = nn, data = 0)
diag(A)                  <- bb
A[cbind(1:(nn-1), 2:nn)] <- cc[-nn]
A[cbind(2:nn, 1:(nn-1))] <- aa[-1]
B <- runif(nn)

# solve as full matrix
solve(A, B)                           

# same,  now using tridiagonal algorithm
as.vector(Solve.tridiag(aa[-1], bb, cc[-nn], B))

# same, now with 3 different right hand sides
B3 <- cbind(B, B*2, B*3)
Solve.tridiag(aa[-1], bb, cc[-nn], B3)

}

\details{
  If the length of the vector \code{dia} is equal to N, then the lengths of 
  \code{diam1} and \code{diap1} should be equal to N-1
}
\seealso{
  \code{\link{Solve.banded}}, the function to solve a banded system of
  linear equations.

  \code{\link{Solve.block}}, the function to solve a block diagonal system of
  linear equations.

  \code{\link{Solve}} the generalised inverse solution,

  \code{\link{solve}} the R default
}
\keyword{array}

