\name{predict.locfit}
\alias{predict.locfit}
\title{
Prediction from a Locfit object.
}
\usage{
\method{predict}{locfit}(object, newdata=NULL, where = "fitp",
          se.fit=FALSE, band="none", what="coef", \dots)
}
\description{
  The \code{\link{locfit}} function computes a local fit at a selected set
of points (as defined by the \code{ev} argument). The \code{predict.locfit}
function is used to interpolate from these points to any other points.
The method is based on cubic hermite polynomial interpolation, using the
estimates and local slopes at each fit point.

  The motivation for this two-step procedure is computational speed.
Depending on the sample size, dimension and fitting procedure, the
local fitting method can be expensive, and it is desirable to keep the
number of points at which the direct fit is computed to a minimum.
The interpolation method used by \code{predict.locfit()} is usually
much faster, and can be computed at larger numbers of points.
}
\arguments{
  \item{object}{Fitted object from \code{\link{locfit}()}.}
  \item{newdata}{Points to predict at. Can be given in several forms:
    vector/matrix; list, data frame.}
  \item{se.fit}{If \code{TRUE}, standard errors are computed along with
    the fitted values.}
  \item{where, what, band}{arguments passed on to
    \code{\link{preplot.locfit}}.} 
  \item{...}{Additional arguments to \code{\link{preplot.locfit}}.}
}
\value{
If \code{se.fit=F}, a numeric vector of predictors.
If \code{se.fit=T}, a list with components \code{fit}, \code{se.fit} and
\code{residual.scale}.
}
\examples{
data(ethanol, package="locfit")
fit <- locfit(NOx ~ E, data=ethanol)
predict(fit,c(0.6,0.8,1.0))
}
\keyword{smooth}
