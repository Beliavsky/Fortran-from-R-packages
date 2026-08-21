\name{kappa0}
\alias{kappa0}
\title{
Critical Values for Simultaneous Confidence Bands.
}
\usage{
kappa0(formula, cov=0.95, ev=lfgrid(20), \dots)
}
\arguments{
\item{formula}{
  Local regression model formula. A \code{"locfit"} object can also
  be provided; in this case the formula and other arguments are extracted
  from this object.
}
\item{cov}{ Coverage Probability for critical values.  }
\item{ev}{ Locfit evaluation structure. Should usually be a grid --
  this specifies the integration rule. }
\item{\dots}{Other arguments to \code{\link{locfit}}. Important arguments
  include \code{flim} and \code{alpha}. }
}
\description{
  The geometric constants for simultaneous confidence bands are computed,
  as described in Sun and Loader (1994) (bias adjustment is not implemented
  here). These are then passed to the \code{\link{crit}} function, which
  computes the critical value for the confidence bands.

  The method requires both the weight diagrams l(x), the
  derivative l'(x) and (in 2 or more dimensions) the second
  derivatives l''(x).
  These are implemented exactly for a constant bandwidth.
  For nearest
  neighbor bandwidths, the computations are approximate and a warning
  is produced.

  The theoretical justification for the bands uses normality of
  the random errors \eqn{e_1,\dots,e_n} in the regression model,
  and in particular the spherical symmetry of the error vector.
  For non-normal distributions, and likelihood models, one relies
  on central limit and related theorems.

  Computation uses the product Simpson's rule to evaluate the
  multidimensional integrals (The domain of integration, and
  hence the region of simultaneous coverage, is determined by
  the \code{flim} argument). Expect the integration to be slow in more
  than one dimension. The \code{mint} argument controls the
  precision.
}
\value{
  A list with components for the critical value, geometric constants,
  e.t.c. Can be passed directly to \code{\link{plot.locfit}} as the
  \code{crit} argument.
}
\seealso{
\code{\link{locfit}}, \code{\link{plot.locfit}},
\code{\link{crit}}, \code{\link{crit<-}}.
}
\examples{
# compute and plot simultaneous confidence bands
data(ethanol)
fit <- locfit(NOx~E,data=ethanol)
crit(fit) <- kappa0(NOx~E,data=ethanol)
plot(fit,crit=crit,band="local")
}
\references{
Sun, J. and Loader, C. (1994). Simultaneous confidence bands for linear
  regression and smoothing. Annals of Statistics 22, 1328-1345.
}
\keyword{smooth}
% Converted by Sd2Rd version 0.2-a5.
