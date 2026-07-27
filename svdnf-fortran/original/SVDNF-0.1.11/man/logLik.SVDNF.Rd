\name{logLik.SVDNF}
\alias{logLik.SVDNF}
\alias{logLik.DNFOptim}
\title{Extract Log-Likelihood for \code{SVDNF} and \code{DNFOptim} Objects}
\usage{
\method{logLik}{SVDNF}(object, ...)
}
\arguments{
 \item{object}{an object of class \code{SVDNF} or \code{DNFOptim}.}
 \item{...}{further arguments passed to or from other methods.}
}
\description{
  Returns the log-likelihood value of the stochastic volatility
  model represented by \code{object} evaluated at the parameters given in the \code{DNF} function.
}
\value{
  The log-likelihood of the stochastic volatility model given by
  \code{object} evaluated at the parameters given to the \code{DNF} function. For \code{DNFOptim} objects, this returns the log-likelihood at the MLE parameter values and the number of free parameters in the model. 
}

\note{It will always be the case \code{df = NULL} for \code{SVDNF} objects as they are evaluations of the DNF algorithm for a fixed set of parameters. However, for \code{DNFOptim} objects, \code{df} will be the number of free parameters in the optimization.}

\examples{
set.seed(1)
# Generate 200 returns from the DuffiePanSingleton model
DuffiePanSingleton_mod <- dynamicsSVM(model = "DuffiePanSingleton") 
DuffiePanSingleton_sim <- modelSim(t = 200, dynamics = DuffiePanSingleton_mod) 

# Run DNF on the data
dnf_filter <- DNF(data = DuffiePanSingleton_sim$returns,
  dynamics = DuffiePanSingleton_mod) 

# Print log-likelihood evaluation.
logLik(dnf_filter)
}
