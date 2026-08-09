\name{varranges}
\alias{varranges}
\title{
  Calculates ranges of inverse variables in a linear inverse problem
}
\description{
  Given the linear constraints

  \deqn{Ex=f}
  \deqn{Gx>=h}
  and a set of "variables" described by the linear equations
  \deqn{Var = EqA.x+EqB}
  finds the minimum and maximum values of the variables
  by successively minimising and maximising each variable equation

}
\usage{
varranges(E=NULL, F=NULL, G=NULL, H=NULL, EqA, EqB=NULL, 
  ispos=FALSE, tol=1e-8, verbose=TRUE, lower=NULL, upper=NULL)
}

\arguments{
  \item{E }{numeric matrix containing the coefficients of the equalities
    \eqn{Ex=F}.
  }
  \item{F }{numeric vector containing the right-hand side of the
    equalities.
  }
  \item{G }{numeric matrix containing the coefficients of the inequalities
    \eqn{Gx>=H}.
  }
  \item{H }{numeric vector containing the right-hand side of the
    inequalities.
  }
  \item{EqA }{numeric matrix containing the coefficients that define the
    variable equations.
  }
  \item{EqB }{numeric vector containing the right-hand side of the variable
    equations.
  }
  \item{ispos }{if \code{TRUE}, it is imposed that unknowns are positive
    quantities.
  }
  \item{tol }{tolerance for equality and inequality constraints.
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
  a 2-column matrix with the minimum and maximum value of each equation
  (variable)
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\examples{
# Ranges in the contribution of food 3+4+5 in the diet of Mink (try ?Minkdiet)

E     <- rbind(Minkdiet$Prey, rep(1, 7))
F     <- c(Minkdiet$Mink, 1)
EqA   <- c(0, 0, 1, 1, 1, 0, 0)                    # sum of food 3,4,5
(isoA <- varranges(E, F, EqA = EqA, ispos = TRUE)) # ranges of part of food 3+4+5

# The same, but explicitly imposing positivity
varranges(E, F, EqA = EqA, G = diag(7), H = rep(0, 7)) 

# The same, but shorter - using lower bound:
varranges(E, F, EqA = EqA, lower=0) 

}
\note{
  uses linear programming function \code{\link[lpSolve]{lp}} from package
  \code{lpSolve}.
}
\references{
Michel Berkelaar and others (2010). lpSolve:
  Interface to Lp_solve v. 5.5 to solve linear/integer
  programs. R package version 5.6.5.
  http://CRAN.R-project.org/package=lpSolve
}
\seealso{
  \code{\link{Minkdiet}}, for a description of the Mink diet example.

  \code{\link{xranges}}, to estimate ranges of inverse unknowns.

  \code{\link{xsample}}, to randomly sample the lsei problem

  \code{\link[lpSolve]{lp}}: linear programming function from package lpSolve.
}

\keyword{algebra}
\keyword{array}
\keyword{optimize}

