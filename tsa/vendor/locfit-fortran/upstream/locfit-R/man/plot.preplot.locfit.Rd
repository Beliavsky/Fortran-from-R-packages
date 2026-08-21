\name{plot.preplot.locfit}
\alias{plot.preplot.locfit}
\title{
Plot a "preplot.locfit" object.
}
\usage{
\method{plot}{preplot.locfit}(x, pv, tv, ...)
}
\arguments{
\item{x}{A \code{preplot.locfit} object, produced by
\code{\link{preplot.locfit}()}.
}
\item{pv, tv, ...}{
Other arguments to \code{plot.locfit.1d}, \code{plot.locfit.2d} or
\code{plot.locfit.3d} as appropriate.
}}
\description{
The \code{\link{plot.locfit}()} function is implemented, roughly, as
a call to \code{\link{preplot.locfit}()}, followed by a call to
\code{plot.locfitpred()}. For most users, there will be little
need to call \code{plot.locfitpred()} directly.
}
\seealso{
\code{\link{locfit}}, \code{\link{plot.locfit}},
\code{\link{preplot.locfit}}, \code{\link{plot.locfit.1d}},
\code{\link{plot.locfit.2d}}, \code{\link{plot.locfit.3d}}.
}
\keyword{smooth}
