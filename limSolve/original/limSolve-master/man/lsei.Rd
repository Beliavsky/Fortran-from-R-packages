\name{lsei}
\alias{lsei}
\title{
  Least Squares with Equalities and Inequalities
}
\description{
  Solves an lsei inverse problem (Least Squares with Equality and Inequality
  Constraints)

  \deqn{\min(||Ax-b||^2)} subject to
  \deqn{Ex=f}
  \deqn{Gx>=h}
 
  Uses either subroutine lsei (FORTRAN) from the LINPACK package, or
  \code{solve.QP} from \R-package \code{quadprog}.

  In case the equality constraints \eqn{Ex=f} cannot be satisfied, a
  generalized inverse solution residual vector length is obtained for \eqn{f-Ex}.

  This is the minimal length possible for \eqn{||f-Ex||^2}.
}

\usage{
lsei (A = NULL, B = NULL, E = NULL, F = NULL, G = NULL, H = NULL,
      Wx = NULL, Wa = NULL, type = 1, tol = sqrt(.Machine$double.eps),
      tolrank = NULL, fulloutput = FALSE, verbose = TRUE, 
      lower = NULL, upper = NULL)
      }

\arguments{
  \item{A }{numeric matrix containing the coefficients of the quadratic
    function to be minimised, \eqn{||Ax-B||^2}; if the columns of \code{A} have
    a names attribute, they will be used to label the output.
  }
  \item{B }{numeric vector containing the right-hand side of the quadratic
    function to be minimised.
  }
  \item{E }{numeric matrix containing the coefficients of the equality
    constraints, \eqn{Ex=F}; if the columns of \code{E} have a names attribute,
    and the columns of \code{A} do not, they will be used to label the output.
  }
  \item{F }{numeric vector containing the right-hand side of the equality
    constraints.
  }
  \item{G }{numeric matrix containing the coefficients of the inequality
    constraints, \eqn{Gx>=H}; if the columns of \code{G} have a names
    attribute, and the columns of \code{A} and \code{E} do not, they will be
    used to label the output.
  }
  \item{H }{numeric vector containing the right-hand side of the inequality
    constraints.
  }
  \item{Wx }{numeric vector with weighting coefficients of unknowns (length
    = number of unknowns).
  }
  \item{Wa }{numeric vector with weighting coefficients of the quadratic
    function (Ax-B) to be minimised (length = number of number of rows of A).
  }
  \item{type }{integer code determining algorithm to use 1=\code{lsei},
    2=\code{solve.QP} from R-package \code{quadprog} (see note).
  }
  \item{tol }{tolerance (for singular value decomposition, equality and
    inequality constraints).
  }
  \item{tolrank }{only used if \code{type} = 1; if not \code{NULL} then
    \code{tolrank} should be a two-valued vector containing the rank
    determination tolerance for the equality constraint equations
    (1st value) and for the reduced least squares equations (2nd value).
  }
  \item{fulloutput }{if \code{TRUE}, also returns the covariance matrix of
    the solution and the rank of the equality constraints - only used
    if \code{type} = 1.
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
  \item{X }{vector containing the solution of the least squares problem.
  }
  \item{residualNorm }{scalar, the sum of absolute values of residuals of
    equalities and violated inequalities.
  }
  \item{solutionNorm }{scalar, the value of the minimised quadratic function
    at the solution, i.e. the value of \eqn{||Ax-b||^2}.
  }
  \item{IsError }{logical, \code{TRUE} if an error occurred.
  }
  \item{type }{the string "lsei", such that how the solution was obtained
    can be traced.
  }
  \item{covar }{covariance matrix of the solution; only returned if
    \code{fulloutput} = \code{TRUE}.
  }
  \item{RankEq }{rank of the equality constraint matrix.; only returned if
    \code{fulloutput} = \code{TRUE}.
  }
  \item{RankApp }{rank of the reduced least squares problem (approximate
    equations); only returned if \code{fulloutput} = \code{TRUE}.
  }
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
# ------------------------------------------------------------------------------
# example 1: polynomial fitting
# ------------------------------------------------------------------------------
x <- 1:5
y <- c(9, 8, 6, 7, 5)
plot(x, y, main = "Polynomial fitting, using lsei", cex = 1.5, 
     pch = 16, ylim = c(4, 10))

# 1-st order
A <- cbind(rep(1, 5), x)
B <- y
cf <- lsei(A, B)$X
abline(coef = cf)

# 2-nd order
A <- cbind(A, x^2)
cf <- lsei(A, B)$X
curve(cf[1] + cf[2]*x + cf[3]*x^2, add = TRUE, lty = 2)

# 3-rd order
A <- cbind(A, x^3)
cf <- lsei(A, B)$X
curve(cf[1] + cf[2]*x + cf[3]*x^2 + cf[4]*x^3, add = TRUE, lty = 3)

# 4-th order
A <- cbind(A, x^4)
cf <- lsei(A, B)$X
curve(cf[1] + cf[2]*x + cf[3]*x^2 + cf[4]*x^3 + cf[5]*x^4, 
      add = TRUE, lty = 4)
legend("bottomleft", c("1st-order", "2nd-order","3rd-order","4th-order"),
       lty = 1:4)

# ------------------------------------------------------------------------------
# example 2: equalities, approximate equalities and inequalities
# ------------------------------------------------------------------------------

A <- matrix(nrow = 4, ncol = 3,
            data = c(3, 1, 2, 0, 2, 0, 0, 1, 1, 0, 2, 0))
B <- c(2, 1, 8, 3)
E <- c(0, 1, 0)
F <- 3
G <- matrix(nrow = 2, ncol = 3, byrow = TRUE,
            data = c(-1, 2, 0, 1, 0, -1))
H <- c(-3, 2)

lsei(E = E, F = F, A = A, B = B, G = G, H = H)

# ------------------------------------------------------------------------------
# example 3: bounds added
# ------------------------------------------------------------------------------
# 
lsei(E = E, F = F, A = A, B = B, G = G, H = H, 
     lower = c(2, NA, 0))

# this cannot be solved
lsei(E = E, F = F, A = A, B = B, G = G, H = H, 
     upper = 2)

# this can
lsei(E = E, F = F, A = A, B = B, G = G, H = H, 
     upper = 3)
}

\references{
  K. H. Haskell and R. J. Hanson, An algorithm for linear least squares problems
  with equality and nonnegativity constraints, Report SAND77-0552,
  Sandia Laboratories, June 1978.

  K. H. Haskell and R. J. Hanson, Selected algorithms for the linearly
  constrained least squares problem - a users guide,
  Report SAND78-1290, Sandia Laboratories,August 1979.
  
  K. H. Haskell and R. J. Hanson, An algorithm for linear least squares
  problems with equality and nonnegativity constraints,
  Mathematical Programming 21 (1981), pp. 98-118.

  R. J. Hanson and K. H. Haskell, Two algorithms for the linearly constrained
  least squares problem, ACM Transactions on Mathematical Software,
  September 1982.

  Berwin A. Turlach R and Andreas Weingessel (2007). quadprog: Functions to
  solve Quadratic Programming Problems. R package version 1.4-11. S original
  by Berwin A. Turlach R port by Andreas Weingessel.
}

\seealso{
  \code{\link{ldei}}, \code{\link{linp}},

  \code{solve.QR} the original function from package \code{quadprog}.
}

\note{
  See comments in the original code for more details; these comments are
  included in the \file{docs} subroutine of the package.
}
\keyword{algebra}
\keyword{array}
\keyword{optimize}

