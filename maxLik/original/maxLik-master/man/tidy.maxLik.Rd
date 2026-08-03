\name{tidy.maxLik}
\alias{tidy.maxLik}
\alias{glance.maxLik}
\title{tidy and glance methods for maxLik objects}
\description{
  These methods return summary information about the estimated model.
  Both require the \pkg{tibble} package to be installed.
}
\usage{
\method{tidy}{maxLik}(x,  ...)
\method{glance}{maxLik}(x, ...)
}
\arguments{
  \item{x}{
    object of class 'maxLik'.
  }
  \item{\ldots}{Not used.}
}
\value{
  For \code{tidy()}, a tibble with columns:
  \describe{
    \item{term}{The name of the estimated parameter (parameters are
      sequentially numbered if names missing).}
    \item{estimate}{The estimated parameter.}
    \item{std.error}{The standard error of the estimate.}
    \item{statistic}{The \eqn{z}{z}-statistic of the estimate.}
    \item{p.value}{The \eqn{p}{p}-value.}
  }
  This is essentially the same table as \code{summary}-method prints,
  just in form of a tibble (data frame).
  
  For \code{glance()}, a one-row tibble with columns:
    \describe{
      \item{df}{The degrees of freedom of the model.}
      \item{logLik}{The log-likelihood of the model.}
      \item{AIC}{Akaike's Information Criterion for the model.}
      \item{nobs}{The number of observations, if this is available, otherwise \code{NA}.}
    }
}
\seealso{
  The functions \code{\link[generics:tidy]{tidy}} and
  \code{\link[generics:glance]{glance}} in package \pkg{generics}, and
    \code{\link[=summary.maxLik]{summary}} to display the
    \dQuote{standard} summary information.
}
\author{David Hugh-Jones}
\examples{
## Example with a single parameter
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t
a <- maxLik(loglik, start=2)
tidy(a)
glance(a)
## Example with a parameter vector
x <- rnorm(100)
loglik <- function(theta) {
   dnorm(x, mean=theta[1], sd=theta[2], log=TRUE)
}
a <- maxLik(loglik, start=c(mu=0, sd=1))
tidy(a)
glance(a)
}
