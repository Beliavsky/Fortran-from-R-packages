\name{nParam}
\alias{nParam}
\alias{nParam.default}
\alias{nParam.lm}
\title{Number of model parameters}
\description{
  This function returns the number of model parameters.  The default
  method returns the component \code{x$param$nParam}.
}
\usage{
nParam(x, free=FALSE, \dots)
\method{nParam}{default}(x, \dots)
\method{nParam}{lm}(x, \dots)
}

\arguments{
  \item{x}{a statistical model}
  \item{free}{logical, whether to report only the free parameters or the
    total number of parameters (default)}
  \item{\dots}{other arguments for methods}
}
\details{
  Free parameters are the parameters with no equality
  restrictions.  Some parameters may be restricted (e.g. sum of two
  probabilities may be restricted to equal unity).   In this case the
  total number of parameters may depend on the
  normalisation.
}
\value{
  Number of parameters in the model
}
\author{Ott Toomet, \email{otoomet@econ.au.dk}}
\seealso{\code{\link{nObs}} for number of observations}
\examples{
# Construct a simple OLS regression:
x1 <- runif(100)
x2 <- runif(100)
y <- 3 + 4*x1 + 5*x2 + rnorm(100)
m <- lm(y~x1+x2)  # estimate it
summary(m)
nParam(m) # you get 3
}
\keyword{methods}
