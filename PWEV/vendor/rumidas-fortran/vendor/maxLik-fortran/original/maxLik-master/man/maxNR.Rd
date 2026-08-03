\name{maxNR}
\alias{maxNR}
\alias{maxBFGSR}
\alias{maxBHHH}
\title{Newton- and Quasi-Newton Maximization}
\description{
  Unconstrained and equality-constrained
  maximization based on the quadratic approximation
  (Newton) method.
  The Newton-Raphson, BFGS (Broyden 1970, Fletcher 1970, Goldfarb 1970,
  Shanno 1970), and BHHH (Berndt, Hall, Hall, Hausman 1974) methods
  are available.
}
\usage{
maxNR(fn, grad = NULL, hess = NULL, start,
      constraints = NULL, finalHessian = TRUE, bhhhHessian=FALSE,
      fixed = NULL, activePar = NULL, control=NULL, ... )
maxBFGSR(fn, grad = NULL, hess = NULL, start,
      constraints = NULL, finalHessian = TRUE,
      fixed = NULL, activePar = NULL, control=NULL, ... )
maxBHHH(fn, grad = NULL, hess = NULL, start, 
      finalHessian = "BHHH", ... )
}
\arguments{
  \item{fn}{the function to be maximized.
    It must have the parameter vector as the first argument and
    it must return either a single number, or a numeric vector (this is
    is summed internally).
    If the BHHH method is used and argument \code{gradient} is not given,
    \code{fn} must return a numeric vector of observation-specific
    log-likelihood values.
    If the parameters are out of range, \code{fn} should
    return \code{NA}.  See details for constant parameters.

    \code{fn} may also return attributes "gradient" and/or "hessian".
    If these attributes are set, the algorithm uses the corresponding
    values as
    gradient and Hessian.
  }
  \item{grad}{gradient of the objective function.
    It must have the parameter vector as the first argument and
    it must return either a gradient vector of the objective function,
    or a matrix, where \emph{columns} correspond to individual parameters.
    The column sums are treated as gradient components.
    If \code{NULL}, finite-difference gradients are computed.
    If BHHH method is used, \code{grad} must return a matrix,
    where rows corresponds to the gradient vectors for individual
    observations and the columns to the individual parameters.
    If \code{fn} returns an object with attribute \code{gradient},
    this argument is ignored.
  }
  \item{hess}{Hessian matrix of the function.
    It must have the parameter vector as the first argument and
    it must return the Hessian matrix of the objective function.
    If missing, finite-difference Hessian, based on \code{gradient},
    is computed.
    Hessian is used by the Newton-Raphson method only, and eventually by
    the other methods if \code{finalHessian} is requested.}
  \item{start}{initial parameter values.  If start values
    are named, those names are also carried over to the results.}
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
    \code{FALSE} (do not calculate), \code{TRUE} (use analytic/finite-difference
    Hessian) or \code{"bhhh"}/\code{"BHHH"} for the information equality
    approach.  The latter approach is only suitable for maximizing
    log-likelihood functions.  It requires the gradient/log-likelihood to
    be supplied by individual observations.
    Note that computing the (actual, not BHHH) final Hessian
    does not carry any extra penalty for the NR method,
    but does for the other methods.}
  \item{bhhhHessian}{logical. Indicating whether to use the information
    equality approximation (Bernd, Hall, Hall, and Hausman, 1974) for
    the Hessian.  This effectively transforms \code{maxNR} into
    \code{maxBHHH} and is mainly designed for internal use.}
  \item{fixed}{parameters to be treated as constants at their
    \code{start} values.  If present, it is treated as an index vector of
    \code{start} parameters.}
  \item{activePar}{this argument is retained for backward compatibility only;
    please use argument \code{fixed} instead.}
  \item{control}{list of control parameters.  The control parameters
    used by these optimizers are
    \describe{
      \item{tol}{\eqn{10^{-8}}{1e-8},
	stopping condition.  Stop if the absolute difference
	between successive iterations is less than \code{tol}.  Return
	\code{code=2}.

	If set to a negative value, the criterion is never fulfilled,
	and hence disabled.
      }
      \item{reltol}{sqrt(.Machine$double.eps), stopping
	condition.  Relative convergence
	tolerance: the algorithm stops if the relative improvement
	between iterations is less than \sQuote{reltol}.  Return code
	8.  Negative value disables condition.
      }
      \item{gradtol}{stopping condition.  Stop if norm of the gradient is
	less than \code{gradtol}.  Return code 1.  Negative value
	disables condition.}
      \item{steptol}{1e-10, stopping/error condition.
	If \code{qac == "stephalving"} and the quadratic
	approximation leads to a worse, instead of a better value, or to
	\code{NA}, the step length
	is halved and a new attempt is made.  If necessary, this procedure is repeated
	until step < \code{steptol}, thereafter code 3 is returned.}
      \item{lambdatol}{\eqn{10^{-6}}{1e-6}, 
	controls whether Hessian is treated as negative
	definite.  If the
	largest of the eigenvalues of the Hessian is larger than
	\code{-lambdatol} (Hessian is not negative definite),
	a suitable diagonal matrix is subtracted from the
	Hessian (quadratic hill-climbing) in order to enforce negative
	definiteness.
      }
      \item{qrtol}{\eqn{10^{-10}}{1e-10},
	QR-decomposition tolerance for the Hessian inversion.
      }
      \item{qac}{"stephalving", Quadratic Approximation Correction.  When the new
	guess is worse than the initial one, the algorithm attemts to correct it:
	"stephalving" decreases the
	step but keeps the direction,
	"marquardt" uses
	\cite{Marquardt (1963)} method by decreasing the step length while also
	moving closer to the pure gradient direction.  It may be faster and
	more robust choice in areas where quadratic approximation
	behaves poorly.  \code{maxNR} and \code{maxBHHH} only.
      }
      \item{marquardt_lambda0}{\eqn{10^{-2}}{1e-2},
	positive numeric, initial correction term for \cite{Marquardt (1963)}
	correction.
      }
      \item{marquardt_lambdaStep}{2, how much the \cite{Marquardt
	  (1963)}
	correction term is
	decreased/increased at each
	successful/unsuccesful step.
	\code{maxNR} and \code{maxBHHH} only.
      }
      \item{marquardt_maxLambda}{\eqn{10^{12}}{1e12},
	maximum allowed \cite{Marquardt (1963)} correction term.  If exceeded, the
	algorithm exits with return code 3.
	\code{maxNR} and \code{maxBHHH} only.
      }
      \item{iterlim}{stopping condition.  Stop if more than \code{iterlim}
	iterations, return \code{code=4}.}
      \item{printLevel}{this argument determines the level of
	printing which is done during the optimization process. The default
	value 0 means that no printing occurs, 1 prints the
	initial and final details, 2 prints all the
	main tracing information for every iteration.  Higher
	values will result in even more output.
      }
    }
  }
  \item{\dots}{further arguments to \code{fn}, \code{grad} and
    \code{hess}.
    Further arguments to \code{maxBHHH} are also passed to
    \code{maxNR}.
    To maintain compatibility with the earlier versions, \dots also passes a
    number of control options (\code{tol}, \code{reltol},
    \code{gradtol}, \code{steptol},
    \code{lambdatol},  \code{qrtol}, \code{iterlim}) to the optimizers.
    }
  }
  \details{
    The idea of the Newton method is to approximate the function at a given
    location by a multidimensional quadratic function, and use the estimated
    maximum as the start value for the next iteration.  Such an
    approximation requires knowledge of both gradient and Hessian, the
    latter of which can be quite costly to compute.  Several methods for
    approximating Hessian exist, including BFGS and BHHH.

The BHHH (information equality) approximation is only valid for
log-likelihood functions.
It requires the score (gradient) values by individual observations and hence
those must be returned 
by individual observations by \code{grad} or \code{fn}.
The Hessian is approximated as the negative of the sum of the outer products
of the gradients of individual observations, or, in the matrix form,
\deqn{
  \mathsf{H}^{BHHH}
  =
  -\frac{1}{N} \sum_{i=1}^N
   \left[
    \frac{\partial \ell(\boldsymbol{\vartheta})}
    {\boldsymbol{\vartheta}}
    \frac{\partial \ell(\boldsymbol{\vartheta})}
    {\boldsymbol{\vartheta}'}
  \right]
  }{
    \code{H = -t(gradient) \%*\% gradient = - crossprod( gradient )}.
    }

The functions \code{maxNR}, \code{maxBFGSR}, and \code{maxBHHH}
can work with constant parameters, useful if a parameter value
converges to the boundary of support, or for testing.  
One way is to put
\code{fixed} to non-NULL, specifying which parameters should be treated as
constants.  The
parameters can also be fixed in runtime (only for \code{maxNR} and \code{maxBHHH}) by
signaling it with the 
\code{fn} return value.  See Henningsen & Toomet (2011) for details.
}

\value{
  object of class "maxim".  Data can be extracted through the following
  methods: 
  \item{\code{\link{maxValue}}}{\code{fn} value at maximum (the last calculated value
    if not converged.)}
  \item{\code{\link[=coef.maxim]{coef}}}{estimated parameter value.}
  \item{gradient}{vector, last calculated gradient value.  Should be
    close to 0 in case of normal convergence.}
  \item{estfun}{matrix of gradients at parameter value \code{estimate}
      evaluated at each observation (only if \code{grad} returns a matrix
      or \code{grad} is not specified and \code{fn} returns a vector).}
  \item{hessian}{Hessian at the maximum (the last calculated value if
    not converged).}
  \item{returnCode}{return code:
    \describe{
    \item{1}{ gradient close to zero (normal convergence).}
    \item{2}{ successive function values within tolerance limit (normal
      convergence).}
    \item{3}{ last step could not find higher value (probably not
      converged).  This is related to line search step getting too
      small, usually because hitting the boundary of the parameter
      space.  It may also be related to attempts to move to a wrong
      direction because of numerical errors.  In some cases it can be
      helped by changing \code{steptol}.}
    \item{4}{ iteration limit exceeded.}
    \item{5}{infinite value.}
    \item{6}{infinite gradient.}
    \item{7}{infinite Hessian.}
    \item{8}{successive function values within relative tolerance
      limit (normal convergence).}
    \item{9}{(BFGS) Hessian approximation cannot be improved because of
      gradient did not change.  May be related to numerical
      approximation problems or wrong analytic gradient.}
    \item{100}{ Initial value out of range.}
    }
  }
  \item{returnMessage}{ a short message, describing the return code.}
  \item{activePar}{logical vector, which parameters are optimized over.
    Contains only \code{TRUE}-s if no parameters are fixed.}
  \item{nIter}{number of iterations.}
  \item{maximType}{character string, type of maximization.}
  \item{maxControl}{the optimization control parameters in the form of a
    \code{\linkS4class{MaxControl}} object.}

  The following components can only be extracted directly (with \code{\$}):
  \item{last.step}{a list describing the last unsuccessful step if
    \code{code=3} with following components:
    \describe{
    \item{theta0}{ previous parameter value}
    \item{f0}{ \code{fn} value at \code{theta0}}
    \item{climb}{ the movement vector to the maximum of the quadratic approximation}
    }
  }
  \item{constraints}{A list, describing the constrained optimization
    (\code{NULL} if unconstrained).  Includes the following components:
    \describe{
      \item{type}{ type of constrained optimization}
      \item{outer.iterations}{ number of iterations in the constraints step}
      \item{barrier.value}{ value of the barrier function}
    }
  }
}

\section{Warning}{
    No attempt is made to ensure that user-provided analytic
    gradient/Hessian is correct.  The users are
    encouraged to use \code{\link{compareDerivatives}} function,
    designed for this purpose.  If analytic gradient/Hessian are wrong,
    the algorithm may not converge, or may converge to a wrong point.

    As the BHHH method
    uses the likelihood-specific information equality,
    it is only suitable for maximizing log-likelihood functions!

    Quasi-Newton methods, including those mentioned above, do not work
    well in non-concave regions.  This is especially the case with the
    implementation in \code{maxBFGSR}.  The user is advised to
    experiment with various tolerance options to achieve convergence.
  }
\references{
   Berndt, E., Hall, B., Hall, R. and Hausman, J. (1974):
   Estimation and Inference in Nonlinear Structural Models,
   \emph{Annals of Social Measurement} \bold{3}, 653--665.

   Broyden, C.G. (1970):
   The Convergence of a Class of Double-rank Minimization Algorithms,
   \emph{Journal of the Institute of Mathematics and Its Applications} \bold{6},
   76--90.

   Fletcher, R. (1970):
   A New Approach to Variable Metric Algorithms,
   \emph{Computer Journal} \bold{13}, 317--322.

   Goldfarb, D. (1970):
   A Family of Variable Metric Updates Derived by Variational Means,
   \emph{Mathematics of Computation} \bold{24}, 23--26.

   Henningsen, A. and Toomet, O. (2011): maxLik: A package for maximum likelihood
   estimation in R \emph{Computational Statistics} \bold{26}, 443--458
   
   Marquardt, D.W., (1963) An Algorithm for Least-Squares Estimation of
   Nonlinear Parameters, \emph{Journal of the Society for Industrial &
   Applied Mathematics} \bold{11}, 2, 431--441

   Shanno, D.F. (1970):
   Conditioning of Quasi-Newton Methods for Function Minimization,
   \emph{Mathematics of Computation} \bold{24}, 647--656.
}

\author{Ott Toomet, Arne Henningsen,
   function \code{maxBFGSR} was originally developed by Yves Croissant
   (and placed in 'mlogit' package)}

\seealso{\code{\link{maxLik}} for a general framework for maximum likelihood
   estimation (MLE);
   \code{\link{maxBHHH}} for maximizations using the Berndt, Hall, Hall,
   Hausman (1974) algorithm (which is a wrapper function to \code{maxNR});
   \code{\link{maxBFGS}} for maximization using the BFGS, Nelder-Mead (NM),
   and Simulated Annealing (SANN) method (based on \code{\link{optim}}),
   also supporting inequality constraints;
   \code{\link{nlm}} for Newton-Raphson optimization; and
  \code{\link{optim}} for different gradient-based optimization
  methods.}

\examples{
## Fit exponential distribution by ML
t <- rexp(100, 2)  # create data with parameter 2
loglik <- function(theta) sum(log(theta) - theta*t)
## Note the log-likelihood and gradient are summed over observations
gradlik <- function(theta) sum(1/theta - t)
hesslik <- function(theta) -100/theta^2
## Estimate with finite-difference gradient and Hessian
a <- maxNR(loglik, start=1, control=list(printLevel=2))
summary(a)
## You would probably prefer 1/mean(t) instead ;-)

## The same example with analytic gradient and Hessian
a <- maxNR(loglik, gradlik, hesslik, start=1)
summary(a)

## BFGS estimation with finite-difference gradient
a <- maxBFGSR( loglik, start=1 )
summary(a)

## For the BHHH method we need likelihood values and gradients
## of individual observations, not the sum of those
loglikInd <- function(theta) log(theta) - theta*t
gradlikInd <- function(theta) 1/theta - t
## Estimate with analytic gradient
a <- maxBHHH(loglikInd, gradlikInd, start=1)
summary(a)

## Example with a vector argument:  Estimate the mean and
## variance of a random normal sample by maximum likelihood
## Note: you might want to use maxLik instead
loglik <- function(param) {
                           # param is a 2-vector of c(mean, sd)
  mu <- param[1]
  sigma <- param[2]
  ll <- -0.5*N*log(2*pi) - N*log(sigma) - sum(0.5*(x - mu)^2/sigma^2)
  ll
}
x <- rnorm(100, 1, 2) # use mean=1, sd=2
N <- length(x)
res <- maxNR(loglik, start=c(0,1)) # use 'wrong' start values
summary(res)

## The previous example with named parameters and a fixed value
resFix <- maxNR(loglik, start=c(mu=0, sigma=1), fixed="sigma")
summary(resFix)  # 'sigma' is exactly 1.000 now.

### Constrained optimization
###
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
res <- maxNR(hatf, start=c(0,0), constraints=list(eqA=A, eqB=B),
             control=list(printLevel=1))
print(summary(res))
}
\keyword{optimize}
