\name{ldp}
\alias{ldp}
\title{
  Least Distance Programming
}

\description{
  Solves the following inverse problem:
  \deqn{\min(\sum {x_i}^2)} subject to
  \deqn{Gx>=h}

  uses least distance programming subroutine ldp (FORTRAN) from Linpack
}

\usage{
ldp(G, H, tol = sqrt(.Machine$double.eps), verbose = TRUE, 
    lower = NULL, upper = NULL)
}

\arguments{
  \item{G }{numeric matrix containing the coefficients of the inequality
    constraints \eqn{Gx>=H}; if the columns of \code{G} have a names attribute,
    they will be used to label the output.
  }
  \item{H }{numeric vector containing the right-hand side of the inequality
    constraints.
  }
  \item{tol }{tolerance (for inequality constraints).
  }
  \item{verbose }{logical to print warnings and messages.
  }
  \item{upper, lower }{vector containing upper and lower bounds 
    on the unknowns. If one value, it is assumed to apply to all unknowns.
    If a vector, it should have a length equal to the number of unknowns; this
    vector can contain NA for unbounded variables. 
    The upper and lower bounds are added to the inequality conditions G*x>=H.
  }
}
\value{
  a list containing:
  \item{X }{vector containing the solution of the least distance problem.
  }
  \item{residualNorm }{scalar, the sum of absolute values of residuals of
    violated inequalities; should be zero or very small if the problem is
    feasible.
  }
  \item{solutionNorm }{scalar, the value of the quadratic function at the
    solution, i.e. the value of \eqn{\sum {w_i*x_i}^2}.
  }
  \item{IsError }{logical, \code{TRUE} if an error occurred.
  }
  \item{type }{the string "ldp", such that how the solution was obtained
    can be traced.
  }
  \item{numiter }{the number of iterations.
  }
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
# parsimonious (simplest) solution
G <- matrix(nrow = 2, ncol = 2, data = c(3, 2, 2, 4))
H <- c(3, 2)

ldp(G, H)

# imposing bounds on the first unknown
ldp(G, H, lower = c(1, NA))
}

\references{
  Lawson C.L.and Hanson R.J. 1974. Solving Least Squares Problems, Prentice-Hall

  Lawson C.L.and Hanson R.J. 1995. Solving Least Squares Problems.
  SIAM classics in applied mathematics, Philadelphia.        (reprint of book)
}

\seealso{
  \code{\link{ldei}}, which includes equalities.
}
\keyword{algebra}
\keyword{array}
\keyword{optimize}

