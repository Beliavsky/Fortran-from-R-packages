\name{lfgrid}
\alias{lfgrid}
\title{
Locfit - grid evaluation structure.
}
\usage{
lfgrid(mg=10, ll, ur)
}
\description{
  \code{lfgrid()} is used to specify evaluation on a grid of points
  for \code{\link{locfit.raw}()}. The structure computes
  a bounding box for the data, and divides that into a grid with
  specified margins.
}
\arguments{
\item{mg}{
Number of grid points along each margin. Can be a single number (which
is applied in each dimension), or a vector specifying a value for
each dimension.
}
\item{ll}{
Lower left limits for the grid. Length should be the number
of dimensions of the data provided to \code{\link{locfit.raw}()}.
}
\item{ur}{
Upper right limits for the grid. By default, \code{ll} and
\code{ur} are generated as the bounding box for the data.
}
}

\examples{
data(ethanol, package="locfit")
plot.eval(locfit(NOx ~ lp(E, C, scale=TRUE), data=ethanol, ev=lfgrid()))
}

\keyword{smooth}
