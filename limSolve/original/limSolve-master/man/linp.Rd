\name{linp}
\alias{linp}
\title{
  Linear Programming.
}

\description{
  Solves a linear programming problem,
  \deqn{\min(\sum {Cost_i.x_i})}{min(sum(Cost_i x_i))} subject to
  \deqn{Ex=f}
  \deqn{Gx>=h}
  \deqn{x_i>=0} (optional)

  This function provides a wrapper around \code{\link[lpSolve]{lp}} (see note)
  from package lpSolve, written to be consistent with the functions
  \code{lsei}, and \code{ldei}.
  
  It allows for the x's to be negative (not standard in lp).
}
\note{
  If the requirement of nonnegativity are relaxed, then strictly speaking the
  problem is not a linear programming problem.

  The function \code{lp} may fail and terminate R for very small problems that
  are repeated frequently...

  Also note that sometimes multiple solutions exist for the same problem.
}

\usage{
linp(E = NULL, F = NULL, G = NULL, H = NULL, Cost,
     ispos = TRUE, int.vec = NULL, verbose = TRUE, 
     lower = NULL, upper = NULL, ...)
}

\arguments{
  \item{E }{numeric matrix containing the coefficients of the equality
    constraints \eqn{Ex=F}; if the columns of \code{E} have a names attribute,
    they will be used to label the output.
  }
  \item{F }{numeric vector containing the right-hand side of the equality
    constraints.
  }
  \item{G }{numeric matrix containing the coefficients of the inequality
    constraints \eqn{Gx>=H}; if the columns of \code{G} have a names attribute,
    and the columns of \code{E} do not, they will be used to label the output.
  }
  \item{H }{numeric vector containing the right-hand side of the inequality
    constraints.
  }
  \item{Cost }{numeric vector containing the coefficients of the cost function;
    if \code{Cost} has a names attribute, and neither the columns of \code{E}
    nor \code{G} have a name, they will be used to label the output.
  }
  \item{ispos}{logical, when \code{TRUE} then the unknowns (x) must be
    positive (this is consistent with the original definition of a
    linear programming problem).
  }
  \item{int.vec }{when not \code{NULL}, a numeric vector giving the indices
    of variables that are required to be an integer. The length of this
    vector will therefore be the number of integer variables.
  }
  \item{verbose }{logical to print warnings and messages.
  }
  \item{upper, lower }{vector containing upper and lower bounds 
    on the unknowns. If one value, it is assumed to apply to all unknowns.
    If a vector, it should have a length equal to the number of unknowns; this
    vector can contain NA for unbounded variables. 
    The upper and lower bounds are added to the inequality conditions G*x>=H.
  }
  \item{... }{extra arguments passed to R-function \code{lp}.
  }
}

\value{
  a list containing:
  \item{X }{vector containing the solution of the linear programming problem.
  }
  \item{residualNorm }{scalar, the sum of absolute values of residuals of
    equalities and violated inequalities.
    Should be very small or zero for a feasible linear programming problem.
  }
  \item{solutionNorm }{scalar, the value of the minimised \code{Cost} function,
    i.e. the value of \eqn{\sum {Cost_i.x_i}}.
  }
  \item{IsError }{logical, \code{TRUE} if an error occurred.
  }
  \item{type }{the string "linp", such that how the solution was obtained can
    be traced.
  }
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\examples{
#-------------------------------------------------------------------------------
# Linear programming problem 1, not feasible
#-------------------------------------------------------------------------------

# maximise x1 + 3*x2 
# subject to
#-x1 -x2    < -3
#-x1 + x2   <-1
# x1 + 2*x2 < 2
# xi > 0

G    <- matrix(nrow = 3, data = c(-1, -1, 1, -1, 1, 2))
H    <- c(3, -1, 2)
Cost <- c(-1, -3)
(L <- linp(E = NULL, F = NULL, Cost = Cost, G = G, H = H))
L$residualNorm

#-------------------------------------------------------------------------------
# Linear programming problem 2, feasible
#-------------------------------------------------------------------------------

# minimise x1 + 8*x2 + 9*x3 + 2*x4 + 7*x5 + 3*x6  
# subject to:
#-x1            + x4 + x5      = 0
#    - x2       - x4      + x6 = 0
# x1 + x2 + x3                 > 1
#           x3       + x5 + x6 < 1
# xi > 0

E <- matrix(nrow = 2, byrow = TRUE, data = c(-1, 0, 0, 1, 1, 0,
                                              0,-1, 0, -1, 0, 1))
F <- c(0, 0)
G <- matrix(nrow = 2, byrow = TRUE, data = c(1, 1, 1, 0, 0, 0,
                                             0, 0, -1, 0, -1, -1))
H    <- c(1, -1)
Cost <- c(1, 8, 9, 2, 7, 3)
(L <- linp(E = E, F = F, Cost = Cost, G = G, H = H))
L$residualNorm

# Including a lower bound:
linp(E = E, F = F, Cost = Cost, G = G, H = H, lower = 0.25)

#-------------------------------------------------------------------------------
# Linear programming problem 3, no positivity
#-------------------------------------------------------------------------------
# minimise x1 + 2x2 -x3 +4 x4
# subject to:
# 3x1 + 2x2 + x3 + x4 = 2
#  x1 +  x2 + x3 + x4 = 2

# 2x1 +  x2 + x3 + x4 >=-1
# -x1 + 3x2 +2x3 + x4 >= 2
# -x1       + x3      >= 1

E <- matrix(ncol = 4, byrow = TRUE,
            data =c(3, 2, 1, 4, 1, 1, 1, 1))
F <- c(2, 2)

G <- matrix(ncol = 4, byrow = TRUE,
            data = c(2, 1, 1, 1, -1, 3, 2, 1, -1, 0, 1, 0))
H <- c(-1, 2, 1)
Cost <- c(1, 2, -1, 4)

linp(E = E, F = F, G = G, H = H, Cost, ispos = FALSE)
}
\seealso{
  \code{\link{ldei}},  \code{\link{lsei}},
  
  \code{\link[lpSolve]{lp}} the original function from package lpSolve
  
  \code{\link{Blending}}, a linear programming problem.
}

\references{
  Michel Berkelaar and others (2007).
  lpSolve: Interface to Lpsolve v. 5.5 to solve linear or integer programs.
  R package version 5.5.8.
}

\keyword{algebra}
\keyword{array}
\keyword{optimize}

