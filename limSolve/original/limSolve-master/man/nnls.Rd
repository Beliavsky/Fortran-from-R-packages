\name{nnls}
\alias{nnls}
\title{
  Nonnegative Least Squares
}
\description{
  Solves the following inverse problem:
  \deqn{\min(||Ax-b||^2)} subject to
  \deqn{x>=0}

  Uses subroutine nnls (FORTRAN) from Linpack
}

\usage{
nnls(A, B, tol = sqrt(.Machine$double.eps), verbose = TRUE)
}

\arguments{
  \item{A }{numeric matrix containing the coefficients of the equality
    constraints \eqn{Ax~=B}; if the columns of \code{A} have a \code{names} 
    attribute, the names will be used to label the output.
  }
  \item{B }{numeric vector containing the right-hand side of the
    equality constraints.
  }
  \item{tol }{tolerance (for singular value decomposition and for the
    "equality" constraints).
  }
  \item{verbose }{logical to print \code{nnls} error messages.
  }
}

\value{
  a list containing:
  \item{X }{vector containing the solution of the nonnegative least
    squares problem.
  }
  \item{residualNorm }{scalar, the sum of absolute values of residuals of
    violated inequalities (i.e. sumof x[<0]); should be zero or very small
    if the problem is feasible.
  }
  \item{solutionNorm }{scalar, the value of the quadratic function at the
    solution, i.e. the value of \eqn{\min(||Ax-b||^2)}.
  }
  \item{IsError }{logical, \code{TRUE} if an error occurred.
  }
  \item{type }{the string "nnls", such that how the solution was obtained
    can be traced.
  }
  \item{numiter }{the number of iterations.
  }
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
A <- matrix(nrow = 2, ncol = 3, data = c(3, 2, 2, 4, 2, 1))
B <- c(-4, 3)
nnls(A, B)
}

\references{
  Lawson C.L.and Hanson R.J. 1974. Solving Least Squares Problems,
  Prentice-Hall

  Lawson C.L.and Hanson R.J. 1995. Solving Least Squares Problems.
  SIAM classics in applied mathematics, Philadelphia.  (reprint of book)
}
\seealso{
  \code{\link{ldei}}, which includes equalities
}

\keyword{algebra}
\keyword{array}
\keyword{optimize}

