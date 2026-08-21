\name{locfit}
\alias{locfit}
\title{
Local Regression, Likelihood and Density Estimation.
}
\usage{
locfit(formula, data=sys.frame(sys.parent()), weights=1, cens=0, base=0,
       subset, geth=FALSE, \dots, lfproc=locfit.raw)
}
\description{
  \code{locfit} is the model formula-based interface to the Locfit
  library for fitting local regression and likelihood models.

  \code{locfit} is implemented as a front-end to \code{\link{locfit.raw}}.
  See that function for options to control smoothing parameters,
  fitting family and other aspects of the fit.
}
\arguments{
\item{formula}{
Model Formula; e.g. \code{y~lp(x)} for a regression model; \code{~lp(x)} for a
density estimation model. Use of \code{lp()} on the RHS is recommended, especially
when non-default smoothing parameters are used.
}
\item{data}{
Data Frame.
}
\item{weights}{
Prior weights (or sample sizes) for individual observations. This is
typically used where observations have unequal variance.
}
\item{cens}{
Censoring indicator. \code{1} (or \code{TRUE}) denotes a censored observation.
\code{0} (or \code{FALSE}) denotes uncensored.
}
\item{base}{
Baseline for local fitting. For local regression models, specifying
a \code{base} is equivalent to using \code{y-base} as the reponse. But \code{base}
also works for local likelihood.
}
\item{subset}{
Subset observations in the data frame.
}
\item{geth}{
Don't use.
}
\item{...}{
Other arguments to \code{\link{locfit.raw}()} (or the \code{lfproc}).
}
\item{lfproc}{
A processing function to compute the local fit. Default is
\code{locfit.raw()}. Other choices include \code{locfit.robust()},
\code{locfit.censor()} and \code{locfit.quasi()}.
}}
\value{
An object with class \code{"locfit"}. A standard set of methods for printing,
ploting, etc. these objects is provided.
}
\seealso{
\code{\link{locfit.raw}}
}
\examples{
# fit and plot a univariate local regression
data(ethanol, package="locfit")
fit <- locfit(NOx ~ E, data=ethanol)
plot(fit, get.data=TRUE)

# a bivariate local regression with smaller smoothing parameter
fit <- locfit(NOx~lp(E,C,nn=0.5,scale=0), data=ethanol)
plot(fit)

# density estimation
data(geyser, package="locfit")
fit <- locfit( ~ lp(geyser, nn=0.1, h=0.8))
plot(fit,get.data=TRUE)
}
\references{
Loader, C. (1999). Local Regression and Likelihood. Springer, New York.
}
\keyword{smooth}
% Converted by Sd2Rd version 0.2-a5.
