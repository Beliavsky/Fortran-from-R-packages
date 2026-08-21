\name{scb}
\alias{scb}
\title{
Simultaneous Confidence Bands
}
\usage{
scb(x, ..., ev = lfgrid(20), simul = TRUE, type = 1)
}
\description{
  \code{scb} is implemented as a front-end to \code{\link{locfit}},
  to compute simultaneous confidence bands using the tube formula
  method and extensions, based on Sun and Loader (1994).
}
\arguments{
\item{x}{A numeric vector or matrix of predictors (as in
  \code{\link{locfit.raw}}), or a model formula (as in
  \code{\link{locfit}}).}
\item{...}{Additional arguments to \code{\link{locfit.raw}}.}
\item{ev}{The evaluation structure to use.  See \code{\link{locfit.raw}}.}
%\item{mg}{
%  The \code{scb()} function evaluates the confidence bands on a grid
%  of points, rather than the default structures used by \code{\link{locfit}}.
%  \code{mg} controls the number of grid points. Default 10.
%}
%\item{flim}{As in  \code{\link{locfit.raw}}, this defaults to the
%  interval (or bounding box, in more than one dimension)  covering
%  the data. The confidence bands are simultaneous over this interval.}
\item{simul}{Should the coverage be simultaneous or pointwise?}
\item{type}{Type of confidence bands. \code{type=0} computes pointwise
  95\% bands. \code{type=1} computes basic simultaneous bands with no
  corrections. \code{type=2,3,4} are the centered and corrected bands
  for parametric regression models listed in Table 3 of Sun, Loader and
  McCormick (2000).}
}
\value{
A list containing the evaluation points, fit, standard deviations and upper
and lower confidence bounds. The class is \code{"scb"}; methods for
printing and ploting are provided.
}
\seealso{
\code{\link{locfit}}, \code{\link{print.scb}}, \code{\link{plot.scb}}.
}
\examples{
# corrected confidence bands for a linear logistic model
data(insect)
fit <- scb(deaths~lp(lconc,deg=1), type=4, w=nins,
           data=insect,family="binomial",kern="parm")
plot(fit)
}
\references{
Sun J. and Loader, C. (1994).
 Simultaneous confidence bands in linear regression and smoothing.
 \emph{The Annals of Statistics} 22, 1328-1345.

Sun, J., Loader, C. and McCormick, W. (2000).
  Confidence bands in generalized linear models.
 \emph{The Annals of Statistics} 28, 429-460.
}
\keyword{smooth}
