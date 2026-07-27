\name{summary.DNFOptim}
\alias{summary.DNFOptim}
\alias{print.summary.DNFOptim}
\title{Summarizing Stochastic Volatility Model Fits from the Discrete Nonlinear Filter}
\usage{
\method{summary}{DNFOptim}(object, confidence, ...)

\method{print}{summary.DNFOptim}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{object}{an object of class \code{DNFOptim} that you want to summary the parameter estimates.}
  \item{x}{an object of class \code{summary.DNFOptim}.}
  \item{confidence}{Confidence level for computing confidence intervals. Should be between 0 and 1. Default is 0.95.}

    \item{digits}{The number of digits to be printed in the \code{print} method for \code{summary.DNFOptim} objects.}
  \item{...}{further arguments passed to or from other methods.}
}
\description{
  Summary method for \code{DNFOptim} objects.
}
\details{
Returns the summary of the output of a \code{DNFOptim} object.
}
\value{
Returns a list with the model used and its estimated parameters.
  \item{model}{The model that was estimated with the \code{DNFOptim} object.}
  \item{coefficients}{Table with the maximum likelihood parameters estimates. If \code{hessian = TRUE} was set in the \code{DNFOptim} function, standard errors and 95\% confidence intervals are given. Then, the table has the following columns:
     (1) \code{Estimate}: the parameter estimate,
     (2) \code{Std Error}: the standard error of the estimate,
      (3) \code{Lower Bound}: the lower bound of the confidence interval, and
     (4) \code{Upper Bound}:The upper bound of the confidence interval.
    }
\item{logLik}{Log-likelihood value at the parameter maximum likelihood estimates and the model's degrees of freedom}
}

\examples{
## For examples see example(DNFOptim)
}
