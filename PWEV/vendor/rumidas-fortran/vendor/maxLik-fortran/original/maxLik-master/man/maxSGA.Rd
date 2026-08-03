\name{maxSGA}
\alias{maxSGA}
\alias{maxAdam}
\title{Stochastic Gradient Ascent}
\description{
  Stochastic Gradient Ascent--based optimizers
}
\usage{
maxSGA(fn = NULL, grad = NULL, hess = NULL, start,
       nObs,
       constraints = NULL, finalHessian = FALSE, 
       fixed = NULL, control=NULL, ... )
maxAdam(fn = NULL, grad = NULL, hess = NULL, start,
        nObs,
        constraints = NULL, finalHessian = FALSE, 
        fixed = NULL, control=NULL, ... )
}
\arguments{
  \item{fn}{the function to be maximized.  As the objective function
    values are not directly used for optimization, this argument is
    optional, given \code{grad} is provided.
    It must have the parameter vector as the first argument, and it must
    have an argument \code{index} to specify the integer index of the selected
    observations.  
    It must return either a single number, or a numeric vector (this is
    is summed internally).
    If the parameters are out of range, \code{fn} should
    return \code{NA}.  See details for constant parameters.

    \code{fn} may also return attributes "gradient" and/or "hessian".
    If these attributes are set, the algorithm uses the corresponding
    values as
    gradient and Hessian.
  }
  \item{grad}{gradient of the objective function.
    It must have the parameter vector as the first argument, and it must
    have an argument \code{index} to specify the integer index of selected
    observations.
    It must return either a gradient vector of the objective function,
    or a matrix, where columns correspond to individual parameters.
    The column sums are treated as gradient components.
    If \code{NULL}, finite-difference gradients are computed.
    If \code{fn} returns an object with attribute \code{gradient},
    this argument is ignored.

    If \code{grad} is not supplied, it is computed by finite-difference
    method using \code{fn}.  However, this is only adviseable for
    small-scale tests, not for any production run.  Obviously, \code{fn}
    must be correctly defined in that case.
  }
  \item{hess}{Hessian matrix of the function.  Mainly for compatibility
    reasons, only used for computing the final Hessian if asked to do
    so by setting \code{finalHessian} to \code{TRUE}. 
    It must have the parameter vector as the first argument and
    it must return the Hessian matrix of the objective function.
    If missing, either finite-difference Hessian, based on
    \code{gradient} or BHHH approach
    is computed if asked to do so.
  }
  \item{start}{initial parameter values.  If these have names, the
    names are also used for results.}
  \item{nObs}{number of observations.  This is used to partition the data
    into individual batches.  The resulting batch
    indices are forwarded to the \code{grad} function through the
    argument \code{index}.}
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
    approach.  The latter approach is only suitable when working with a
    log-likelihood function, and it requires the gradient/log-likelihood to
    be supplied by individual observations.

    Hessian matrix is not often used for optimization problems where one
    applies SGA, but even if one is not interested in standard errors,
    it may provide useful information about the model performance.  If
    computed by finite-difference method, the Hessian computation may be
    very slow.
  }
  \item{fixed}{parameters to be treated as constants at their
    \code{start} values.  If present, it is treated as an index vector of
    \code{start} parameters.}
  \item{control}{list of control parameters.  The ones
    used by these optimizers are
    \describe{
      \item{SGA_momentum}{0, numeric momentum parameter for SGA.  Must lie
	in interval \eqn{[0,1]}{[0,1]}.  See details.
      }
      Adam-specific parameters
      \item{Adam_momentum1}{0.9, numeric in interval \eqn{(0,1)}{(0,1)}, the first moment momentum}
      \item{Adam_momentum2}{0.999, numeric in interval \eqn{(0,1)}{(0,1)}, the second moment momentum}
      General stochastic gradient parameters:
      \item{SG_learningRate}{step size the SGA algorithm takes in the
	gradient direction.  If 1, the step equals to the gradient value.  A
	good value is often 0.01--0.3}
      \item{SG_batchSize}{SGA batch size, an integer between 1 and
	\code{nObs}.
	If \code{NULL} (default), the full batch gradient is computed.
      }
      \item{SG_clip}{\code{NULL}, gradient clipping threshold.  The
	algorithm ensures that \eqn{||g(\theta)||_2^2 \le
	  \kappa}{norm(gradient)^2 <= kappa} where \eqn{\kappa}{kappa} is
	the \code{SG_clip} value.  If the
	actual norm of the gradient exceeds (square root of)
	\eqn{\kappa}{kappa},
	the gradient will be scaled back accordingly while
	preserving its direction.  \code{NULL} means no clipping.
      }
      Stopping conditions:
      \item{gradtol}{stopping condition.  Stop if norm of the gradient is
	less than \code{gradtol}.  Default 0, i.e. do not use this
	condition.  This condition is useful if the
	objective is to drive full batch gradient to zero on training data.
	It is not a good objective in case of the stochastic
	gradient, and if the objective is to optimize the objective on
	validation data.
      }
      \item{SG_patience}{\code{NULL}, or integer.  Stopping condition:
	the algorithm counts how many times
	the objective function has been worse than its best value so
	far, and if this exceeds \code{SG_patience}, the algorithm stops.
      }
      \item{SG_patienceStep}{1L, integer.  After how many epochs to check
	the patience value.  \code{1} means to check at each epoch, and hence to compute the
	objective function.  This may be undesirable if the objective
	function is costly to compute.
      }
      \item{iterlim}{stopping condition.  Stop if more than \code{iterlim}
	epochs, return \code{code=4}.
	Epoch is a set of iterations that cycles through all
	observations.  In case of full batch, iterations and epochs are
	equivalent.  If \code{iterlim = 0}, does not do any learning and
	returns the initial values unchanged.
      }
      \item{printLevel}{this argument determines the level of
	printing which is done during the optimization process. The default
	value 0 means that no printing occurs, 1 prints the
	initial and final details, 2 prints all the
	main tracing information for every epoch.  Higher
	values will result in even more output.
      }
      \item{storeParameters}{logical, whether to store and return the
	parameter
	values at each epoch.  If \code{TRUE}, the stored values
	can be retrieved with \code{\link{storedParameters}}-method.  The
	parameters are stored as a matrix with rows corresponding to the
	epochs and columns to the parameter components.  There are
	\code{iterlim} + 1 rows, where the first one corresponds to the
	initial parameters.

	Default \code{FALSE}.
      }
      \item{storeValues}{logical, whether to store and return the objective
	function values at each epoch.  If \code{TRUE}, the stored values
	can be retrieved with \code{\link{storedValues}}-method.  There are
	\code{iterlim} + 1 values, where the first one corresponds to
	the value at the
	initial parameters.
	
	Default \code{FALSE}.
      }
    }
    See \code{\link{maxControl}} for more information.
  }
  \item{\dots}{further arguments to \code{fn}, \code{grad} and
    \code{hess}.
    To maintain compatibility with the earlier versions, \dots also
    passes certain control options to the optimizers.
  }
}

\details{
  Gradient Ascent (GA) is a optimization method where the algorithm
  repeatedly takes small steps in the gradient's direction, the
  parameter vector \eqn{\theta}{theta} is updated as \eqn{\theta
    \leftarrow theta + \mathrm{learning rate}\cdot \nabla
    f(\theta)}{theta <- learning rate * gradient f(theta)}.
  In case of Stochastic GA (SGA), the gradient is not computed on the
  full set of observations but on a small subset, \emph{batch},
  potentially a single observation only.  In certain circumstances
  this converges much faster
  than when using all observation (see
  \cite{Bottou et al, 2018}).
  
  If \code{SGA_momentum} is positive, the SGA algorithm updates the parameters
  \eqn{\theta}{theta} in two steps.  First, the momentum is used to update
  the \dQuote{velocity} \eqn{v}{v} as
  \eqn{v \leftarrow \mathrm{momentum}\cdot v + \mathrm{learning
      rate}\cdot \nabla f(\theta)}{v <- momentum*v + learning
    rate* gradient f(theta)}, and thereafter the parameter
  \eqn{\theta}{theta} is updates as
  \eqn{\theta \leftarrow \theta + v}{theta <- theta + v}.  Initial
  velocity is set to 0.
  
  The Adam algorithm is more complex and uses first and second moments
  of stochastic gradients to automatically adjust the learning rate.
  See \cite{Goodfellow et al, 2016, page 301}.
    
  The function \code{fn} is not directly used for optimization, only
  for printing or as a stopping condition.  In this sense
  it is up to the user to decide what the function
  returns, if anything.  For instance, it may be useful for \code{fn} to compute the
  objective function on either full training data, or on validation data,
  and just ignore the \code{index} argument.  The latter is useful if
  using \emph{patience}-based stopping.
  However, one may also
  choose to select the observations determined by the index to
  compute the objective function on the current data batch.
  % Does it support contraints?
}

\value{
  object of class "maxim".  Data can be extracted through the following
  methods: 
  \item{\code{\link{maxValue}}}{\code{fn} value at maximum (the last calculated value
    if not converged.)}
  \item{\code{\link{coef}}}{estimated parameter value.}
  \item{\code{\link{gradient}}}{vector, last calculated gradient value.  Should be
    close to 0 in case of normal convergence.}
  \item{estfun}{matrix of gradients at parameter value \code{estimate}
      evaluated at each observation (only if \code{grad} returns a matrix
      or \code{grad} is not specified and \code{fn} returns a vector).}
  \item{\code{\link{hessian}}}{Hessian at the maximum (the last calculated value if
    not converged).}
  \item{\code{\link{storedValues}}}{return values stored at each epoch}
  \item{\code{\link{storedParameters}}}{return parameters stored at each epoch}
  \item{\code{\link{returnCode}}}{
    a numeric code that describes the convergence or error.
  }
  \item{\code{\link{returnMessage}}}{a short message, describing the return code.}
  \item{\code{\link{activePar}}}{logical vector, which parameters are optimized over.
    Contains only \code{TRUE}-s if no parameters are fixed.}
  \item{\code{\link{nIter}}}{number of iterations.}
  \item{\code{\link{maximType}}}{character string, type of maximization.}
  \item{\code{\link{maxControl}}}{the optimization control parameters in the form of a
    \code{\linkS4class{MaxControl}} object.}
}

\references{
    Bottou, L.; Curtis, F. & Nocedal, J.:
    Optimization Methods for
    Large-Scale Machine Learning \emph{SIAM Review}, 2018, \bold{60},
    223--311.

    Goodfellow, I.; Bengio, Y.; Courville, A. (2016): Deep Learning,
    \emph{MIT Press}

    Henningsen, A. and Toomet, O. (2011): maxLik: A package for maximum likelihood
    estimation in R \emph{Computational Statistics} \bold{26}, 443--458
}

\author{Ott Toomet, Arne Henningsen}

\seealso{
  A good starting point to learn about the usage of stochastic gradient
  ascent in \pkg{maxLik} package is the vignette \dQuote{Stochastic
    Gradient Ascent in maxLik}.

  The other related functions are
   \code{\link{maxNR}} for Newton-Raphson, a popular Hessian-based maximization;
   \code{\link{maxBFGS}} for maximization using the BFGS, Nelder-Mead (NM),
   and Simulated Annealing (SANN) method (based on \code{\link[stats]{optim}}),
   also supporting inequality constraints;
   \code{\link{maxLik}} for a general framework for maximum likelihood
   estimation (MLE);
   \code{\link{optim}} for different gradient-based optimization
   methods.
 }

 \examples{
## estimate the exponential distribution parameter by ML
set.seed(1)
t <- rexp(100, 2)
loglik <- function(theta, index) sum(log(theta) - theta*t[index])
## Note the log-likelihood and gradient are summed over observations
gradlik <- function(theta, index) sum(1/theta - t[index])
## Estimate with full-batch
a <- maxSGA(loglik, gradlik, start=1, control=list(iterlim=1000,
            SG_batchSize=10), nObs=100)
            # note that loglik is not really needed, and is not used
            # here, unless more print verbosity is asked
summary(a)
##
## demonstrate the usage of index, and using
## fn for computing the objective function on validation data.
## Create a linear model where variables are very unequally scaled
##
## OLS loglik function: compute the function value on validation data only
loglik <- function(beta, index) {
   e <- yValid - XValid \%*\% beta
   -crossprod(e)/length(y)
}
## OLS gradient: compute it on training data only
## Use 'index' to select the subset corresponding to the minibatch
gradlik <- function(beta, index)  {
   e <- yTrain[index] - XTrain[index,,drop=FALSE] \%*\% beta
   g <- t(-2*t(XTrain[index,,drop=FALSE]) \%*\% e)
   -g/length(index)
}
N <- 1000
## two random variables: one with scale 1, the other with 100
X <- cbind(rnorm(N), rnorm(N, sd=100))
beta <- c(1, 1)  # true parameter values
y <- X \%*\% beta + rnorm(N, sd=0.2)
## training-validation split
iTrain <- sample(N, 0.8*N)
XTrain <- X[iTrain,,drop=FALSE]
XValid <- X[-iTrain,,drop=FALSE]
yTrain <- y[iTrain]
yValid <- y[-iTrain]
##
## do this without momentum: learning rate must stay small for the gradient not to explode
cat("  No momentum:\n")
a <- maxSGA(loglik, gradlik, start=c(10,10),
           control=list(printLevel=1, iterlim=50,
                        SG_batchSize=30, SG_learningRate=0.0001, SGA_momentum=0
                        ), nObs=length(yTrain))
print(summary(a))  # the first component is off, the second one is close to the true value
## do with momentum 0.99
cat("  Momentum 0.99:\n")
a <- maxSGA(loglik, gradlik, start=c(10,10),
           control=list(printLevel=1, iterlim=50,
                        SG_batchSize=30, SG_learningRate=0.0001, SGA_momentum=0.99
                        # no momentum
                        ), nObs=length(yTrain))
print(summary(a))  # close to true value
}
\keyword{optimize}
