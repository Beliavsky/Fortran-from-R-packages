\name{Solve.banded}
\alias{Solve.banded}
\title{
  Solution of a banded system of linear equations
}
\description{
  Solves the linear system of equations \deqn{Ax = B} by Gaussion elimination

  where \code{A} has to be square, and \emph{banded}, i.e. with the only nonzero
  elements in bands near the diagonal.
 
  The matrix \code{A} is either inputted as a full square matrix or as the non-zero
  bands.

  uses lapack subroutine dgbsv (FORTRAN)
}

\usage{
Solve.banded(abd, nup, nlow, B = rep(0, times = ncol(abd)),
  full = (nrow(abd) == ncol(abd)))}

\arguments{
  \item{abd }{either a matrix containing the (nonzero) bands, rotated row-wise
    (anti-clockwise) only, or a full square matrix.
  }
  \item{nup }{number of nonzero bands above the diagonal; ignored if \code{full}
    matrix is inputted.
  }
  \item{nlow }{number of nonzero bands below the diagonal; ignored if \code{full}
    matrix is inputted.
  }
  \item{B }{Right-hand side of the equations, a vector with length = number
   of rows of \code{A}, or a matrix with number of rows = number of rows 
   of \code{A}.
  }
  \item{full }{if \code{TRUE}: full matrix is inputted,
    if \code{FALSE}: banded matrix is input.
  }
}

\value{
  matrix with the solution, \code{X}, of the banded system of equations A X =B,
  the number of columns of this matrix = number of columns of \code{B}.
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\note{A similar function but that requires a totally different input can now 
  also be found in the \code{Matrix} package
}

\examples{
# 1. Generate a banded matrix of random numbers, full format
nup  <- 2                         # nr nonzero bands above diagonal
ndwn <- 3                         # nr nonzero bands below diagonal
nn   <- 10                        # nr rows and columns of A
A <- matrix(nrow = nn, ncol = nn, data = runif(1 : (nn*nn)))
A [row(A) < col(A) - nup | row(A) > col(A) + ndwn] <- 0 
diag(A) <- 1                      # 1 on diagonal is easily recognised 

# right hand side
B <- runif(nrow(A))                

# solve it, using the default solver and banded (inputting full matrix)
Full  <- solve(A, B)
Band1 <- Solve.banded(A, nup, ndwn, B)

# 2. create banded form of matrix A
Aext <- rbind(matrix(ncol = ncol(A), nrow = nup, 0),
              A, 
              matrix(ncol = ncol(A), nrow = ndwn, 0))
              
abd  <- matrix(nrow = nup + ndwn + 1, ncol = nn,
               data = Aext[col(Aext) <= row(Aext) & 
                           col(Aext) >= row(Aext) - ndwn - nup])

# print both to screen
A
abd

# solve problem with banded version
Band2 <- Solve.banded(abd, nup, ndwn, B)

# compare 3 methods of solution
cbind(Full, Band1, Band2)

# same, now with 3 different right hand sides
B3 <- cbind(B, B*2, B*3)
Solve.banded(abd, nup, ndwn, B3)

}

\references{
  J.J. Dongarra, J.R. Bunch, C.B. Moler, G.W. Stewart,
    LINPACK Users' Guide, SIAM, 1979.
}

\details{
  If the input matrix \code{abd} is square, it is assumed that the full,
  square A is inputted, unless \code{full} is set to \code{FALSE}.
  
	If \code{abd} is not square, then the number of columns denote the
  number of unknowns, while the number of rows equals the nonzero bands,
  i.e. nup+nlow+1
}

\seealso{
  \code{\link{Solve.tridiag}} to solve a tridiagonal system of linear equations.

  \code{\link{Solve}} the generalised inverse solution,

  \code{\link{solve}} the R default
}
\keyword{array}

