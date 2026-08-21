\name{preplot.locfit}
\alias{preplot.locfit}
\title{
Prediction from a Locfit object.
}
\usage{
\method{preplot}{locfit}(object, newdata=NULL, where, tr=NULL, what="coef",
  band="none", get.data=FALSE, f3d=FALSE, \dots)
}
\arguments{
  \item{object}{Fitted object from \code{\link{locfit}()}.}
  \item{newdata}{Points to predict at. Can be given in several forms:
    vector/matrix; list, data frame.}
  \item{where}{An alternative to \code{newdata}. Choices include
    \code{"grid"} for the grid \code{\link{lfmarg}(object)};
    \code{"data"} for the original data points and \code{"fitp"} for the
    direct fitting points (ie. no interpolation).}
  \item{tr}{Transformation for likelihood models. Default is the inverse
    of the link function.}
  \item{what}{What to compute predicted values of. The default,
    \code{what="coef"}, works with the fitted curve itself. Other
    choices include \code{"nlx"} for the length of the weight diagram;
    \code{"infl"} for the influence function; \code{"band"} for the
    bandwidth; \code{"degr"} for the local polynomial degree;
    \code{"lik"} for the maximized local likelihood; \code{"rdf"} for
    the local residual degrees of freedom and \code{"vari"} for the
    variance function. The interpolation algorithm for some of these
    quantities is questionable.}
  \item{band}{Compute standard errors for the fit and include confidence
    bands on the returned object. Default is \code{"none"}. Other
    choices include \code{"global"} for bands using a global variance
    estimate; \code{"local"} for bands using a local variance estimate
    and \code{"pred"} for prediction bands (at present, using a global
    variance estimate).  To obtain the global variance estimate for a
    fit, use \code{\link{rv}}.  This can be changed with
    \code{\link{rv<-}}. Confidence bands, by default, are 95\%, based on
    normal approximations and neglecting bias.  To change the critical
    value or confidence level, or to obtain simultaneous instead of
    pointwise confidence, the critical value stored on the fit must be
    changed. See the \code{\link{kappa0}} and \code{\link{crit}}
    functions.}
  \item{get.data}{If \code{TRUE}, the original data is attached to the
    returned object, and added to the plot.}
  \item{f3d}{If \code{TRUE}, sets a flag that forces ploting using the
    trellis style.  Not available in R.}
  \item{\dots}{arguments passed to and from other methods.}
}
\value{
An object with class \code{"preplot.locfit"}, containing the predicted
values and additional information used to construct the plot.
}
\description{
  \code{preplot.locfit} can be called directly, although it is more usual
  to call \code{\link{plot.locfit}} or \code{\link{predict.locfit}}.
  The advantage of \code{preplot.locfit} is in S-Plus 5, where arithmetic
  and transformations can be performed on the \code{"preplot.locfit"}
  object.

  \code{plot(preplot(fit))} is essentially synonymous with \code{plot(fit)}.
}
\seealso{
\code{\link{locfit}}, \code{\link{predict.locfit}}, \code{\link{plot.locfit}}.
}
\keyword{smooth}
% Converted by Sd2Rd version 0.2-a5.
