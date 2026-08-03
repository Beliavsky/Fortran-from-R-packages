\name{maxLik}
\alias{maxLik}
\alias{print.maxLik}

\title{Maximum likelihood estimation}
\description{
  This is the main interface for the \pkg{maxLik} package, and the
  function that performs Maximum
  Likelihood estimation.  It is a wrapper for different optimizers
  returning an object
  of class "maxLik".  Corresponding methods handle the
  likelihood-specific properties of the estimates,  including standard
  errors. 
}
\usage{
maxLik(logLik, grad = NULL, hess = NULL, start, method,
constraints=NULL, ...)
}
\arguments{
  \item{logLik}{log-likelihood function.  Must have the parameter vector
    as the first argument.  Must return either a single log-likelihood
    value, or a numeric vector where each component is log-likelihood
    of the corresponding individual observation.}
  \item{grad}{gradient of log-likelihood.  Must have the parameter
    vector as the first argument.  Must return either a single gradient
    vector with length equal to the number of parameters, or a matrix
    where each row is the gradient vector of the corresponding individual
    observation.  If \code{NULL}, numeric gradient will be used.}
  \item{hess}{hessian of log-likelihood.  Must have the parameter
    vector as the first argument.  Must return a square matrix.  If
    \code{NULL}, numeric Hessian will be used.}
  \item{start}{numeric vector, initial value of parameters.  If it has
    names, these will also be used for naming the results.}
  \item{method}{maximisation method, currently either 
   "NR" (for Newton-Raphson),
   "BFGS" (for Broyden-Fletcher-Goldfarb-Shanno), 
   "BFGSR" (for the BFGS algorithm implemented in \R),
   "BHHH" (for Berndt-Hall-Hall-Hausman), 
   "SANN" (for Simulated ANNealing), 
   "CG" (for Conjugate Gradients), 
   or "NM" (for Nelder-Mead).  
   Lower-case letters (such as "nr" for Newton-Raphson) are allowed.
   The default method is "NR" for unconstrained problems, and "NM" or
   "BFGS" for constrained problems, depending on if the \code{grad}
   argument was provided.  "BHHH" is a good alternative given the
   likelihood is returned observation-wise (see \code{\link{maxBHHH}}).

   Note that stochastic gradient ascent (SGA) is currently not supported
   as this method seems to be rarely used for maximum likelihood estimation.
 }
  \item{constraints}{either \code{NULL} for unconstrained maximization
    or a list, specifying the constraints.  See \code{\link{maxBFGS}}.
  }
  \item{\dots}{further arguments, such as \code{control},
    \code{iterlim}, or \code{tol},
    are passed to the selected maximisation routine,
    i.e. \code{\link{maxNR}}, \code{\link{maxBFGS}}, \code{\link{maxBFGSR}},
    \code{\link{maxBHHH}}, \code{\link{maxSANN}}, \code{\link{maxCG}},
    or \code{\link{maxNM}}
    (depending on argument \code{method}).  Arguments not used by the
    optimizers are forwarded to \code{logLik}, \code{grad} and
    \code{hess}.
   }
}
\details{
  \code{maxLik} supports constrained optimization in the sense that
  constraints are passed further to the underlying optimization
  routines, and suitable default method is selected.  However, no
  attempt is made to correct the resulting variance-covariance matrix.
  Hence the inference may be wrong.  A corresponding warning is issued
  by the summary method.
}
\value{
  object of class 'maxLik' which inherits from class 'maxim'.
  Useful methods
  include
  \itemize{
    \item \code{\link[=AIC.maxLik]{AIC}}: estimated parameter value
    \item \code{\link[=coef.maxLik]{coef}}: estimated parameter value
    \item \code{\link[=logLik.maxLik]{logLik}}: log-likelihood value
    \item \code{\link{nIter}}: number of iterations
    \item \code{\link[=stdEr.maxLik]{stdEr}}: standard errors
    \item \code{\link[=summary.maxLik]{summary}}: print summary table
    with estimates, standard errors, p, and z-values.
    \item \code{\link[=vcov.maxLik]{vcov}}: variance-covariance matrix
  }
}
\section{Warning}{The constrained maximum likelihood estimation should
  be considered experimental.  In particular, the variance-covariance
  matrix is not corrected for constrained parameter space.
  }
\author{Ott Toomet, Arne Henningsen}
\seealso{\code{\link{maxNR}}, \code{\link{nlm}} and \code{\link{optim}}
  for different non-linear optimisation routines, see
  \code{\link{maxBFGS}} for the constrained maximization examples.}
\examples{
## Estimate the parameter of exponential distribution
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t
gradlik <- function(theta) 1/theta - t
hesslik <- function(theta) -100/theta^2
## Estimate with numeric gradient and hessian
a <- maxLik(loglik, start=1, control=list(printLevel=2))
summary( a )
##
## Estimate with analytic gradient and hessian.
## require much smaller tolerance
## setting 'tol=0' or negative essentially disables this stopping criterion
a <- maxLik(loglik, gradlik, hesslik, start=1,
            control=list(tol=-1, reltol=1e-12, gradtol=1e-12))
summary( a )
##
## Next, we give an example with vector argument:
## fit normal distribution by estimating mean and standard deviation
## by maximum likelihood
##
loglik <- function(param) {
                           # param: vector of 2, c(mean, standard deviation)
   mu <- param[1]
   sigma <- param[2]
   ll <- -0.5*N*log(2*pi) - N*log(sigma) - sum(0.5*(x - mu)^2/sigma^2)
                           # can use dnorm(x, mu, sigma, log=TRUE) instead
   ll
}
x <- rnorm(100, 1, 2) # use mean=1, stdd=2
N <- length(x)
res <- maxLik(loglik, start=c(0,1)) # use 'wrong' start values
summary(res)
##
## Same example, but now with named parameters and a fixed value
##
resFix <- maxLik(loglik, start=c(mu=0, sigma=1), fixed="sigma")
summary(resFix)  # 'sigma' is exactly 1.000 now.
}
\keyword{optimize}
