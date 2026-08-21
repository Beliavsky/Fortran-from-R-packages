\name{uobyqa}
\alias{uobyqa}
\encoding{UTF-8}
\title{An R interface to the uobyqa implementation of Powell}
\description{
  The purpose of \code{uobyqa} is to minimize a function of many variables 
  by a trust region method that forms quadratic models by interpolation.
 
}
\usage{
uobyqa(par, fn, control = list(), \dots)
}
\arguments{
 \item{par}{A numeric vector of starting estimates. }
 \item{fn}{A function that returns the value of the objective at the
   supplied set of parameters \code{par} using auxiliary data in \dots.
   The first argument of \code{fn} must be \code{par}. }
 \item{control}{
    An optional list of control settings.  See the details section for
    the names of the settable control values and their effect.
 }
 \item{\dots}{Further arguments to be passed to \code{fn}.}
}
\details{
  Functions \code{fn} must return a numeric value.

  The \code{control} argument is a list.  Possible named values in the
  list and their defaults are:
  \describe{
    \item{rhobeg}{
      \code{rhobeg} and \code{rhoend} must be set to the initial and final
      values of a trust region radius, so both must be positive with
      \code{0 < rhoend < rhobeg}. Typically \code{rhobeg} should be about
      one tenth of the greatest expected change to a variable.
    }
    \item{rhoend}{
      The smallest value of the trust region radius that is allowed. If
      not defined, then 1e-6 times the value set for \code{rhobeg} will be
      used.
    }
      \item{iprint}{
      The value of \code{iprint} should be set to an integer value in
      \code{0, 1, 2, 3, ...},
      which controls the amount of printing. Specifically, there is no
      output if \code{iprint=0} and there is output only at the start
      and the return if
      \code{iprint=1}. Otherwise, each new value of \code{rho} is printed,
      with the best vector of variables so far and the corresponding value
      of the objective function. Further, each new value of the objective
      function with its variables are output if \code{iprint=3}.
      If \code{iprint > 3}, the objective
      function value and corresponding variables are output every \code{iprint}
      evaluations.
      Default value is \code{0}.
    }
    \item{maxfun}{
      The maximum allowed number of function evaluations. If this is
      exceeded, the method will terminate.
    } 
  }

  Powell's Fortran code has been slightly modified (thanks to Doug Bates for 
  help on this) to avoid use of PRINT statements. Output is now via calls to
  C routines set up to work with the routines BOBYQA, NEWUOA and UOBYQA.
  }
\value{
  A list with components:
  \item{par}{The best set of parameters found.}
  \item{fval}{The value of the objective at the best set of parameters found.}
  \item{feval}{The number of function evaluations used.}
  \item{ierr}{An integer error code.  A value of zero indicates
    success. Other values (consistent with BOBYQA values) are
    \describe{
       \item{1}{maximum number of function evaluations exceeded}
       \item{3}{a trust region step failed to reduce q (Consult Powell for explanation.)}
    }
  } 
  \item{msg}{A message describing the outcome of UOBYQA}
}
\references{
     M. J. D. Powell, "The uobyqa software for unconstrained
     optimization without derivatives",
     in \bold{Large-Scale Nonlinear Optimization},
     Series: Nonconvex Optimization and 
     Its Applications , Vol. 83, Di Pillo, Gianni; Roma, Massimo (Eds.)
     2006, New York: Springer US.

     M. J. D. Powell,
     "Developments of uobyqa for minimization without derivatives",
     IMA Journal of  Numerical Analysis, 2008; 28: 649-664.

     
     Description was taken from comments in the Fortran code of
     M. J. D. Powell on which \pkg{minqa} is based.
}
\seealso{\code{\link{optim}}, \code{\link{nlminb}}}
\examples{
fr <- function(x) {   ## Rosenbrock Banana function
    100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
}
(x3 <- uobyqa(c(1, 2), fr))
## => optimum at c(1, 1) with fval = 0
# check the error exits
# too many iterations
x3e<-uobyqa(c(1, 2), fr, control = list(maxfun=50))
str(x3e)


# To add if we can find them -- examples of ierr = 3.


}
\keyword{nonlinear}
\keyword{optimize}

