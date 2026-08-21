\name{fitted.locfit}
\alias{fitted.locfit}
\title{
Fitted values for a `"locfit"' object.
}
\usage{
\method{fitted}{locfit}(object, data=NULL, what="coef", cv=FALSE,
studentize=FALSE, type="fit", tr, ...)
}
\description{Evaluates the fitted values (i.e. evaluates the surface
  at the original data points) for a Locfit object. This function works
  by reconstructing the model matrix from the original formula, and
  predicting at those points. The function may be fooled; for example,
  if the original data frame has changed since the fit, or if the
  model formula includes calls to random number generators.}
\arguments{
\item{object}{
\code{"locfit"} object.
}
\item{data}{
The data frame for the original fit. Usually, this shouldn't be needed,
especially when the function is called directly. It may be needed
when called inside another function.
}
\item{what}{
What to compute fitted values of. The default, \code{what="coef"}, works
  with the fitted curve itself. Other choices include \code{"nlx"} for the
  length of the weight diagram; \code{"infl"} for the influence function;
  \code{"band"} for the bandwidth; \code{"degr"} for the local polynomial
  degree; \code{"lik"} for the maximized local likelihood; \code{"rdf"}
  for the local residual degrees of freedom and \code{"vari"} for the
  variance function. The interpolation algorithm for some of these quantities
  is questionable.
}
\item{cv}{
If \code{TRUE}, leave-one-out cross validated fitted values are approximated.
Won't make much sense, unless \code{what="coef"}.
}
\item{studentize}{
If \code{TRUE}, residuals are studentized.
}
\item{type}{
  Type of fit or residuals to compute. The default is \code{"fit"} for
  \code{fitted.locfit}, and \code{"dev"} for \code{\link{residuals.locfit}}.
  Other choices include \code{"pear"} for Pearson residuals; \code{"raw"}
  for raw residuals, \code{"ldot"} for likelihood derivative;
  \code{"d2"} for the deviance residual squared; \code{lddot} for the
  likelihood second derivative. Generally, \code{type} should only be
  used when \code{what="coef"}.
}
\item{tr}{
Back transformation for likelihood models.
}
\item{...}{arguments passed to and from methods.}
}
\value{
A numeric vector of the fitted values.
}
\seealso{
  \code{locfit},
  \code{predict.locfit},
  \code{residuals.locfit}
}
\keyword{smooth}
% Converted by Sd2Rd version 0.2-a5.
