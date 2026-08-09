\name{xranges}
\alias{xranges}
\title{
  Calculates ranges of the unknowns of a linear inverse problem
}
\description{
  Given the linear constraints
  \deqn{Ex=f}
  \deqn{Gx>=h}
  finds the minimum and maximum values of all elements of vector \eqn{x}

  This is done by successively minimising and maximising each \code{x}, 
  using linear programming. 
}

\usage{
xranges(E = NULL, F = NULL, G = NULL, H = NULL, 
        ispos = FALSE, tol = 1e-8, central = FALSE, full=FALSE,
        verbose = TRUE, lower = NULL, upper = NULL)}

\arguments{
  \item{E }{numeric matrix containing the coefficients of the equalities
    \eqn{Ex=F}.
  }
  \item{F }{numeric vector containing the right-hand side of the equalities.
  }
  \item{G }{numeric matrix containing the coefficients of the inequalities
    \eqn{Gx>=H}.
  }
  \item{H }{numeric vector containing the right-hand side of the inequalities.
  }
  \item{ispos }{if \code{TRUE}, it is imposed that unknowns are positive
    quantities. Used when generating \code{x0} by linear programming, if x0 is not given.
  }
  \item{tol }{tolerance for equality and inequality constraints.
  }
  \item{central }{if TRUE, the mean value of all range solutions is also
    outputted.
  }
  \item{full }{if \code{TRUE}, all range solutions are also outputted.
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
  a matrix with at least two columns:

  column 1 and 2: the minimum and maximum value of each \code{x}

  if \code{central} is TRUE: column 3 = the central value

  if \code{full} is TRUE: next columns contain all valid range solutions
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
# Estimate the ranges in the Diet Composition of Mink
E <- rbind(Minkdiet$Prey, rep(1, 7))
F <- c(Minkdiet$Mink, 1)
(DC <- xranges(E, F, ispos = TRUE))

# The same, but explicitly imposing positivity using G and H
xranges(E, F, G = diag(7), H = rep(0, 7))

# and using lower bound
xranges(E, F, lower = 0, verbose = FALSE)

}
\details{
  The ranges are estimated by successively minimising and maximising each
  unknown, and using linear programming (based on function \code{lp} from 
  R-package \code{lpSolve}.

  By default linear programming assumes that all unknowns are positive.
  If all unknowns are indeed to be positive, then it will generally be faster
  to set \code{ispos} equal to TRUE
  If \code{ispos} is \code{FALSE}, then a system double the size of the 
  original system must be solved.

  \code{xranges} outputs only the minimum and maximum value of each flow unless:

  \code{full} is \code{TRUE}. In this case, all the results of the successive
  minimisation and maximisation will be outputted, i.e. for each linear
  programming application, not just the value of the unknown being optimised
  but also the corresponding values of the other unknowns will be outputted.

  If \code{central} is \code{TRUE}, then the mean of all the results of the
  linear programming will be outputted.
  This may be a good starting value for \code{\link{xsample}}

  Note: the columns corresponding to the \code{central} value and the
  \code{full} results are valid solutions of the equations \eqn{Ex=F}
  and \eqn{Gx>=H}. This is not the case for the first two columns (with
  the minimal and maximal values).
}
\references{
Michel Berkelaar and others (2010). lpSolve:
  Interface to Lp_solve v. 5.5 to solve linear/integer
  programs. R package version 5.6.5.
  http://CRAN.R-project.org/package=lpSolve
}
\seealso{
  \code{\link{Minkdiet}}, for a description of the Mink diet example.

  \code{\link{varranges}}, for range estimation of variables,

  \code{\link{xsample}}, to randomly sample the lsei problem

 \code{\link[lpSolve]{lp}}: linear programming from package lpSolve
}

\keyword{algebra}
\keyword{array}
\keyword{optimize}

