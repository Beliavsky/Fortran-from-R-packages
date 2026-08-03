\name{sumt}
\Rdversion{1.1}
\alias{sumt}
\title{
  Equality-constrained optimization
}
\description{
  Sequentially Unconstrained Maximization Technique (SUMT) based
  optimization for linear equality constraints.

  This implementation is primarily intended to be called from other
  maximization routines, such as \code{\link{maxNR}}.
}
\usage{
sumt(fn, grad=NULL, hess=NULL,
start,
maxRoutine, constraints,
SUMTTol = sqrt(.Machine$double.eps),
SUMTPenaltyTol = sqrt(.Machine$double.eps),
SUMTQ = 10,
SUMTRho0 = NULL,
printLevel=print.level, print.level = 0, SUMTMaxIter = 100, ...)
}
\arguments{
  \item{fn}{
    function of a (single) vector parameter.  The function may have more
    arguments (passed by \dots), but those are not treated as the
    parameter. 
  }
  \item{grad}{
    gradient function of \code{fn}.  NULL if missing
  }
  \item{hess}{
    function, Hessian of the \code{fn}.  NULL if missing
  }
  \item{start}{
    numeric, initial value of the parameter
}
  \item{maxRoutine}{
    maximization algorithm, such as \code{\link{maxNR}}
}
  \item{constraints}{list, information for constrained maximization.
    Currently two components are supported: \code{eqA} and \code{eqB}
    for linear equality constraints: \eqn{A \beta + B = 0}{A \%*\% beta
      + B = 0}.  The user must ensure that the matrices \code{A} and
    \code{B} are conformable.}
  \item{SUMTTol}{
    stopping condition.  If the estimates at successive outer
    iterations are close enough, i.e. maximum of the absolute value over
    the component difference is smaller than SUMTTol, the algorithm
    stops.

    Note this does not necessarily mean that the constraints are
    satisfied.
    If the penalty function is too \dQuote{weak}, SUMT may repeatedly find
    the same optimum.  In that case a warning is issued.  The user may
    set SUMTTol to a lower value, e.g. to zero.
  }
  \item{SUMTPenaltyTol}{
    stopping condition.  If the barrier value (also called penalty)
    \eqn{(A \beta + B)'(A \beta + B)}{t(A \%*\% beta
      + B) \%*\% (A \%*\% beta + B)} is less than
    \code{SUMTTol}, the algorithm stops
  }
  \item{SUMTQ}{
    a double greater than one, controlling the growth of the \code{rho}
    as described in Details. Defaults to 10.
  }
  \item{SUMTRho0}{
    Initial value for \code{rho}.  If not specified, a (possibly)
    suitable value is
    selected.  See Details.

    One should consider supplying \code{SUMTRho0} in case where the
    unconstrained problem does not have a maximum, or the maximum is too
    far from the constrained value.  Otherwise the authomatically selected
    value may not lead to convergence.
    }
  \item{printLevel}{
    Integer, debugging information.  Larger number prints more details.
  }
  \item{print.level}{same as \sQuote{printLevel}, for backward
    compatibility} 
  \item{SUMTMaxIter}{
    Maximum SUMT iterations
  }
  \item{\dots}{
    Other arguments to \code{maxRoutine} and \code{fn}.
  }
}
\details{
  The Sequential Unconstrained Minimization Technique is a heuristic
  for constrained optimization.  To minimize a function \eqn{f}{f}
  subject to
  constraints, it uses a non-negative penalty function \eqn{P}{P},
  such that \eqn{P(x)}{P(x)} is zero
     iff \eqn{x}{x}
     satisfies the constraints.  One iteratively minimizes
     \eqn{f(x) + \varrho_k P(x)}{f(x) + rho_k P(x)}, where the
     \eqn{\varrho}{rho}
     values are increased according to the rule
     \eqn{\varrho_{k+1} = q \varrho_k}{rho_{k+1} = q rho_k} for some
     constant \eqn{q > 1}{q > 1}, until convergence is
     achieved in the sense that the barrier value
  \eqn{P(x)'P(x)}{P(x)'P(x)} is close to zero.  Note that there is no
     guarantee that the global constrained optimum is
     found.  Standard practice recommends to use the best solution
     found in \dQuote{sufficiently many} replications.

     Any of
     the maximization algorithms in the \pkg{maxLik}, such as
     \code{\link{maxNR}}, can be used for the unconstrained step.

     Analytic gradient and hessian
     are used if provided.
}
\value{
  Object of class 'maxim'.  In addition, a component
  \item{constraints}{A list, describing the constrained optimization.
    Includes the following components:
    \describe{
      \item{type}{type of constrained optimization}
      \item{barrier.value}{value of the penalty function at maximum}
      \item{code}{code for the stopping condition}
      \item{message}{a short message, describing the stopping condition}
      \item{outer.iterations}{number of iterations in the SUMT step}
    }
  }
}
\section{Note}{
  In case of equality constraints, it
  may be more efficient to enclose the function
  in a wrapper function.  The wrapper calculates full set of parameters
  based on a smaller set of parameters, and the
  constraints.
}
\author{
  Ott Toomet, Arne Henningsen
}
\seealso{
  \code{\link[clue]{sumt}} in package \pkg{clue}.
}
\examples{
## We maximize exp(-x^2 - y^2) where x+y = 1
hatf <- function(theta) {
   x <- theta[1]
   y <- theta[2]
   exp(-(x^2 + y^2))
   ## Note: you may prefer exp(- theta \%*\% theta) instead
}
## use constraints: x + y = 1
A <- matrix(c(1, 1), 1, 2)
B <- -1
res <- sumt(hatf, start=c(0,0), maxRoutine=maxNR,
            constraints=list(eqA=A, eqB=B))
print(summary(res))
}
\keyword{optimize}
