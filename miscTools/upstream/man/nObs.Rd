\name{nObs}
\alias{nObs}
\alias{nObs.default}
\alias{nObs.lm}
\title{Return number of observations for statistical models}
\description{
  Returns number of observations for statistical models.  The default
  method assumes presence of a component \code{param$nObs} in \code{x}.
}
\usage{
nObs(x, \dots)
\method{nObs}{default}(x, \dots)
\method{nObs}{lm}(x, \dots)
}
\arguments{
  \item{x}{a statistical model, such as created by \code{\link{lm}}}
  \item{\dots}{further arguments for methods}
}
\details{
  This is a generic function.  The default method returns the component
  \code{x$param$nObs}.  The \code{lm}-method is based on
  qr-decomposition, in the same way as the does \code{\link{summary.lm}}.
}
\value{
  numeric, number of observations
}
\author{Ott Toomet, \email{otoomet@econ.au.dk}}

\examples{
# Construct a simple OLS regression:
x1 <- runif(100)
x2 <- runif(100)
y <- 3 + 4*x1 + 5*x2 + rnorm(100)
m <- lm(y~x1+x2)  # estimate it
nObs(m)
}
\keyword{methods}
