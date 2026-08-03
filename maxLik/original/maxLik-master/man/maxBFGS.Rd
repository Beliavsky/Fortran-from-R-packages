\name{maxBFGS}
\alias{maxBFGS}
\alias{maxCG}
\alias{maxSANN}
\alias{maxNM}

\title{BFGS, conjugate gradient, SANN and Nelder-Mead Maximization}
\description{
  These functions are wrappers for \code{\link{optim}}, adding
  constrained optimization and fixed parameters.
}
\usage{
maxBFGS(fn, grad=NULL, hess=NULL, start, fixed=NULL,
   control=NULL,
   constraints=NULL,
   finalHessian=TRUE,
   parscale=rep(1, length=length(start)),
   ... )

maxCG(fn, grad=NULL, hess=NULL, start, fixed=NULL,
   control=NULL,
   constraints=NULL,
   finalHessian=TRUE,
   parscale=rep(1, length=length(start)), ...)

maxSANN(fn, grad=NULL, hess=NULL, start, fixed=NULL,
   control=NULL,
   constraints=NULL,
   finalHessian=TRUE,
   parscale=rep(1, length=length(start)),
   ... )

maxNM(fn, grad=NULL, hess=NULL, start, fixed=NULL,
   control=NULL,
   constraints=NULL,
   finalHessian=TRUE,
   parscale=rep(1, length=length(start)),
   ...)
}

\arguments{
  \item{fn}{function to be maximised.  Must have the parameter vector as
    the first argument.  In order to use numeric gradient
    and BHHH method, \code{fn} must return a vector of
    observation-specific likelihood values.  Those are summed internally where
    necessary.  If the parameters are out of range, \code{fn} should
    return \code{NA}.  See details for constant parameters.}
  \item{grad}{gradient of \code{fn}.  Must have the parameter vector as
    the first argument.  If \code{NULL}, numeric
    gradient is used (\code{maxNM} and \code{maxSANN} do not use
    gradient).  
    Gradient may return
    a matrix, where columns correspond to the parameters and rows to the
    observations (useful for maxBHHH).  The columns are summed internally.}
  \item{hess}{Hessian of \code{fn}.  Not used by any of these methods, included for
    compatibility with \code{\link{maxNR}}.}
  \item{start}{initial values for the parameters.  If start values
    are named, those names are also carried over to the results.}
  \item{fixed}{parameters to be treated as constants at their
    \code{start} values.  If present, it is treated as an index vector of
    \code{start} parameters.}
  \item{control}{list of control parameters or a \sQuote{MaxControl}
    object.  If it is a list, the default values are used for the
    parameters that are left unspecified by the user.
    These functions accept the following parameters:
    \describe{
      \item{reltol}{sqrt(.Machine$double.eps), stopping
	condition.  Relative convergence
	tolerance: the algorithm stops if the relative improvement
	between iterations is less than \sQuote{reltol}.  Note: for
	compatibility reason \sQuote{tol} is equivalent to
	\sQuote{reltol} for optim-based optimizers.
      }
      \item{iterlim}{integer, maximum number of iterations.
	Default values are 200 for \sQuote{BFGS}, 500
	(\sQuote{CG} and \sQuote{NM}), and 10000 (\sQuote{SANN}).
	Note that
	\sQuote{iteration} may mean different things for different
	optimizers.
      }
      \item{printLevel}{integer, larger number prints more working
	information.  Default 0, no information.
      }
      \item{nm_alpha}{1, Nelder-Mead simplex method reflection
	coefficient (see Nelder & Mead, 1965)
      }
      \item{nm_beta}{0.5, Nelder-Mead contraction coefficient}
      \item{nm_gamma}{2, Nelder-Mead expansion coefficient}
      % SANN
      \item{sann_cand}{\code{NULL} or a function for \code{"SANN"} algorithm
	to generate a new candidate point;
	if \code{NULL}, Gaussian Markov kernel is used
	(see argument \code{gr} of \code{\link{optim}}).}
      \item{sann_temp}{10, starting temperature
	for the \dQuote{SANN} cooling schedule.  See \code{\link{optim}}.}
      \item{sann_tmax}{10, number of function evaluations at each temperature for
	the \dQuote{SANN} optimizer.  See \code{\link{optim}}.}
      \item{sann_randomSeed}{123, integer to seed random numbers to
	ensure replicability of \dQuote{SANN} optimization and preserve
	\code{R} random numbers.  Use
	options like \code{sann_randomSeed=Sys.time()} or
	\code{sann_randomSeed=sample(100,1)} if you want stochastic
	results.
	}
      }
    }
  \item{constraints}{either \code{NULL} for unconstrained optimization
    or a list with two components.  The components may be either
    \code{eqA} and \code{eqB} for equality-constrained optimization
    \eqn{A \theta + B = 0}{A \%*\% theta + B = 0}; or \code{ineqA} and
    \code{ineqB} for inequality constraints \eqn{A \theta + B > 0}{A
      \%*\% theta + B > 0}.  More
       than one
       row in \code{ineqA} and \code{ineqB} corresponds to more than
       one linear constraint, in that case all these must be zero
       (equality) or positive (inequality constraints).
       The equality-constrained problem is forwarded
    to \code{\link{sumt}}, the inequality-constrained case to
    \code{\link{constrOptim2}}.
  }
  \item{finalHessian}{how (and if) to calculate the final Hessian.  Either
  \code{FALSE} (not calculate), \code{TRUE} (use analytic/numeric
  Hessian) or \code{"bhhh"}/\code{"BHHH"} for information equality
  approach.  The latter approach is only suitable for maximizing
  log-likelihood function.  It requires the gradient/log-likelihood to
  be supplied by individual observations, see \code{\link{maxBHHH}} for
  details. 
 }
  \item{parscale}{A vector of scaling values for the parameters.
          Optimization is performed on 'par/parscale' and these should
          be comparable in the sense that a unit change in any element
          produces about a unit change in the scaled value. (see
    \code{\link{optim}})}
  \item{\dots}{further arguments for \code{fn} and \code{grad}.}
}

\details{
  In order to provide a consistent interface, all these functions also
  accept arguments that other optimizers use.  For instance,
  \code{maxNM} accepts the \sQuote{grad} argument despite being a
  gradient-less method.
  
   The \sQuote{state} (or \sQuote{seed}) of R's random number generator
   is saved at the beginning of the \code{maxSANN} function
   and restored at the end of this function
   so this function does \emph{not} affect the generation of random numbers
   although the random seed is set to argument \code{random.seed}
   and the \sQuote{SANN} algorithm uses random numbers.
}

\value{
  object of class "maxim".  Data can be extracted through the following
  functions: 
  \item{maxValue}{\code{fn} value at maximum (the last calculated value
    if not converged.)}
  \item{coef}{estimated parameter value.}
  \item{gradient}{vector, last calculated gradient value.  Should be
    close to 0 in case of normal convergence.}
  \item{estfun}{matrix of gradients at parameter value \code{estimate}
      evaluated at each observation (only if \code{grad} returns a matrix
      or \code{grad} is not specified and \code{fn} returns a vector).}
  \item{hessian}{Hessian at the maximum (the last calculated value if
    not converged).}
  \item{returnCode}{integer. Success code, 0 is success (see
    \code{\link{optim}}).}
  \item{returnMessage}{ a short message, describing the return code.}
  \item{activePar}{logical vector, which parameters are optimized over.
    Contains only \code{TRUE}-s if no parameters are fixed.}
  \item{nIter}{number of iterations.  Two-element integer vector giving the number of
    calls to \code{fn} and \code{gr}, respectively.
    This excludes those calls needed to
    compute the Hessian, if requested, and any calls to \code{fn} to compute a
    finite-difference approximation to the gradient.}
  \item{maximType}{character string, type of maximization.}
  \item{maxControl}{the optimization control parameters in the form of a
    \code{\linkS4class{MaxControl}} object.}

  The following components can only be extracted directly (with \code{\$}):
  \item{constraints}{A list, describing the constrained optimization
    (\code{NULL} if unconstrained).  Includes the following components:
    \describe{
      \item{type}{type of constrained optimization}
      \item{outer.iterations}{number of iterations in the constraints step}
      \item{barrier.value}{value of the barrier function}
    }
  }
}
\author{Ott Toomet, Arne Henningsen}

\seealso{\code{\link{optim}}, \code{\link{nlm}}, \code{\link{maxNR}},
  \code{\link{maxBHHH}}, \code{\link{maxBFGSR}} for a
  \code{\link{maxNR}}-based BFGS implementation.}

\references{
  Nelder, J. A. & Mead, R. A, Simplex Method for Function
  Minimization, The Computer Journal, 1965, 7, 308-313
}

\examples{
# Maximum Likelihood estimation of Poissonian distribution
n <- rpois(100, 3)
loglik <- function(l) n*log(l) - l - lfactorial(n)
# we use numeric gradient
summary(maxBFGS(loglik, start=1))
# you would probably prefer mean(n) instead of that ;-)
# Note also that maxLik is better suited for Maximum Likelihood
###
### Now an example of constrained optimization
###
f <- function(theta) {
  x <- theta[1]
  y <- theta[2]
  exp(-(x^2 + y^2))
  ## you may want to use exp(- theta \%*\% theta) instead
}
## use constraints: x + y >= 1
A <- matrix(c(1, 1), 1, 2)
B <- -1
res <- maxNM(f, start=c(1,1), constraints=list(ineqA=A, ineqB=B),
control=list(printLevel=1))
print(summary(res))
}
\keyword{optimize}
