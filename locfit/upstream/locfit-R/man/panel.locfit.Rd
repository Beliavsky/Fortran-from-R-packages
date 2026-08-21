\name{panel.locfit}
\alias{panel.locfit}
\title{
Locfit panel function
}
\usage{
panel.locfit(x, y, subscripts, z, rot.mat, distance, shade,
                 light.source, xlim, ylim, zlim, xlim.scaled,
                 ylim.scaled, zlim.scaled, region, col, lty, lwd,
                 alpha, col.groups, polynum, drape, at, xlab, ylab,
                 zlab, xlab.default, ylab.default, zlab.default,
                 aspect, panel.aspect, scales.3d, contour, labels,
                 ...)
}
\description{
This panel function can be used to add locfit fits to plots generated
by Lattice.

Currently it works with \code{xyplot} for 1-d fits
and crudely with \code{wireframe} for 2-d fits.
}
\arguments{
  \item{x, y, subscripts, z}{usual arguments to a \code{panel} function}
  \item{rot.mat, distance, shade,
                 light.source, xlim, ylim, zlim, xlim.scaled,
                 ylim.scaled, zlim.scaled, region, col, lty, lwd,
                 alpha, col.groups, polynum, drape, at, xlab, ylab,
                 zlab, xlab.default, ylab.default, zlab.default,
                 aspect, panel.aspect, scales.3d, contour, labels}{
                 further arguments passed on to underlying plotting functions}
  \item{...}{Most Locfit arguments can be passed through \code{xyplot}.}
}
\examples{
\dontrun{
# a simple multi-panel display
data(ethanol, package="locfit")
xyplot(NOx ~ E | C, data=ethanol, panel=panel.locfit)

# The second example uses some Locfit optional arguments.
# Note we can pass the alpha (bandwidth) and family arguments directly to
# xyplot. The cens argument must be given in full; not as a data frame variable.
# The resulting plot does not (yet) distinguish the censored points, but
# the fit will correctly apply censoring.
data(border, package="locfit")
xyplot(runs ~ day, data=border, panel=panel.locfit, family="poisson",
 alpha=0.3, cens=border$no)
}
}
\seealso{
  \code{\link{locfit}}, \code{\link{plot.locfit.3d}}, \code{xyplot}.
}

\keyword{internal}
