\name{lp}
\alias{lp}
\title{
Local Polynomial Model Term
}
\usage{
lp(..., nn, h, adpen, deg, acri, scale, style)
}
\description{
  \code{lp} is a local polynomial model term for Locfit models.
  Usually, it will be the only term on the RHS of the model formula.

  Smoothing parameters should be provided as arguments to \code{lp()},
  rather than to \code{\link{locfit}()}.
}
\arguments{
\item{...}{Predictor variables for the local regression model.
}
\item{nn}{
Nearest neighbor component of the smoothing parameter.
Default value is 0.7, unless either \code{h} or \code{adpen} are
provided, in which case the default is 0.
}
\item{h}{
The constant component of the smoothing parameter. Default: 0.
}
\item{adpen}{Penalty parameter for adaptive fitting.}
\item{deg}{Degree of polynomial to use.}
\item{acri}{Criterion for adaptive bandwidth selection.}
\item{style}{Style for special terms (\code{\link{left}},
  \code{\link{ang}} e.t.c.). Do not try to set this directly;
  call \code{\link{locfit}} instead. }
\item{scale}{
A scale to apply to each variable. This is especially important for
multivariate fitting, where variables may be measured in
non-comparable units. It is also used to specify the frequency
for \code{\link{ang}} terms. If \code{scale=F} (the default) no scaling
is performed. If \code{scale=T}, marginal standard deviations are used.
Alternatively, a numeric vector can provide scales for the
individual variables.
}
}

\seealso{
\code{\link{locfit}},
\code{\link{locfit.raw}}
}
\examples{
data(ethanol, package="locfit")
# fit with 50% nearest neighbor bandwidth.
fit <- locfit(NOx~lp(E,nn=0.5),data=ethanol)
# bivariate fit.
fit <- locfit(NOx~lp(E,C,scale=TRUE),data=ethanol)

# density estimation
data(geyser, package="locfit")
fit <- locfit.raw(lp(geyser,nn=0.1,h=0.8))
}
\keyword{models}
