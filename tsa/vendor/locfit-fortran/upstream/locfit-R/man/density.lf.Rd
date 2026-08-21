\name{density.lf}
\alias{density.lf}
\title{
Density estimation using Locfit
}
\usage{
density.lf(x, n = 50, window = "gaussian", width, from, to,
  cut = if(iwindow == 4.) 0.75 else 0.5,
  ev = lfgrid(mg = n, ll = from, ur = to),
  deg = 0, family = "density", link = "ident", ...)
}
\description{
  This function provides an interface to Locfit, in the syntax of
  (a now old version of) the S-Plus \code{density} function. This can reproduce
  \code{density} results, but allows additional
  \code{\link{locfit.raw}} arguments, such as the degree of fit, to be given.

  It also works in double precision, whereas \code{density} only works
  in single precision.
}
\arguments{
\item{x}{numeric vector of observations whose density is to be estimated.}
\item{n}{number of evaluation points.
  Equivalent to the \code{\link{locfit.raw} mg} argument.}
\item{window}{Window type to use for estimation.
  Equivalent to the \code{\link{locfit.raw} kern} argument. This includes
  all the \code{density} windows except \code{cosine}.}
\item{width}{Window width. Following \code{density}, this is the full
  width; not the half-width usually used by Locfit and
  many other smoothers.}
\item{from}{Lower limit for estimation domain.}
\item{to}{Upper limit for estimation domain.}
\item{cut}{Controls default expansion of the domain.}
\item{ev}{Locfit evaluation structure -- default \code{\link{lfgrid}()}.}
\item{deg}{Fitting degree -- default 0 for kernel estimation.}
\item{family}{Fitting family -- default is \code{"density"}.}
\item{link}{Link function -- default is the \code{"identity"}.}
\item{...}{Additional arguments to \code{\link{locfit.raw}}, with standard
  defaults.}
}
\value{
  A list with components \code{x} (evaluation points) and \code{y}
  (estimated density).
}
\seealso{
\code{density},
\code{\link{locfit}},
\code{\link{locfit.raw}}
}
\examples{
data(geyser)
density.lf(geyser, window="tria")
# the same result with density, except less precision.
density(geyser, window="tria")
}
\keyword{smooth}
