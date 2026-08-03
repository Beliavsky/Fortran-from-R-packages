\name{confint.maxLik}
\alias{confint.maxLik}
\alias{confint}
\title{confint method for maxLik objects}
\description{
  Wald confidence intervals for Maximum Likelihood Estimates
}
\usage{
\method{confint}{maxLik}(object, parm, level=0.95,  ...)
}
\arguments{
  \item{object}{
    object of class \dQuote{maxLik} returned by \code{\link{maxLik}} function
  }
  \item{parm}{the name of parameters to compute the confidence
    intervals.  If omitted, confidence intervals for all parameters are
    computed.}
  \item{level}{the level of confidence interval
  }
  \item{\dots}{additional arguments to be passed to the other methods
  }
}
\value{
  A matrix of lower and upper confidence interval limits (in the first
  and second column respectively).  The matrix rows are labeled by the
  parameter names (if any) and columns by the corresponding distribution
  quantiles. 
}
\seealso{
  \code{\link[stats]{confint}} for the generic \code{confint}
  function, 
  \code{\link[=stdEr.maxLik]{stdEr}} for computing standard errors and
  \code{\link[=summary.maxLik]{summary}} for summary output that includes statistical
  significance information.
}
\author{Luca Scrucca}
\examples{
## compute MLE parameters of normal random sample
x <- rnorm(100)
loglik <- function(theta) {
   dnorm(x, mean=theta[1], sd=theta[2], log=TRUE)
}
m <- maxLik(loglik, start=c(mu=0, sd=1))
summary(m)
confint(m)
confint(m, "mu", level=0.1)
}
