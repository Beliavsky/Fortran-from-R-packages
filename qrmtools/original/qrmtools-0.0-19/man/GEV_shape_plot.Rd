\name{GEV_shape_plot}
\alias{GEV_shape_plot}
\title{Fitted GEV Shape as a Function of the Threshold}
\description{
  Fit GEVs to block maxima and plot the fitted GPD shape as
  a function of the block size.
}
\usage{
GEV_shape_plot(x, blocksize = tail(pretty(seq_len(length(x)/20), n = 64), -1),
               estimate.cov = TRUE, conf.level = 0.95,
               CI.col = adjustcolor(1, alpha.f = 0.2),
               lines.args = list(), xlim = NULL, ylim = NULL,
               xlab = "Block size",  ylab = NULL,
               xlab2 = "Number of blocks", plot = TRUE, ...)
}
\arguments{
  \item{x}{\code{\link{vector}} of \code{\link{numeric}} data.}
  \item{blocksize}{\code{\link{numeric}} vector of block sizes for which
    to fit a GEV to the block maxima.}
  \item{estimate.cov}{\code{\link{logical}} indicating whether
    confidence intervals are to be computed.}
  \item{conf.level}{confidence level of the confidence intervals if
    \code{estimate.cov}.}
  \item{CI.col}{color of the pointwise asymptotic confidence intervals
    (CIs); if \code{NA}, no CIs are shown.}
  \item{lines.args}{\code{\link{list}} of arguments passed to
    the underlying \code{\link{lines}()} for drawing the shape
    parameter as a function of the block size.}
  \item{xlim, ylim, xlab, ylab}{see \code{\link{plot}()}.}
  \item{xlab2}{label of the secondary x-axis.}
  \item{plot}{\code{\link{logical}} indicating whether a plot is produced.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{plot}()}.}
}
\value{
  Invisibly returns a \code{\link{list}} containing the block sizes
  considered, the corresponding block maxima and the fitted GEV
  distribution objects as returned by the underlying
  \code{\link{fit_GEV_MLE}()}.
}
\details{
  Such plots can be used in the block maxima method for determining
  the optimal block size (as the smallest after which the plot
  is (roughly) stable).
}
\author{Marius Hofert}
\examples{
set.seed(271)
X <- rPar(5e4, shape = 4)
GEV_shape_plot(X)
abline(h = 1/4, lty = 3) # theoretical xi = 1/shape for Pareto
}
\keyword{distribution}
\keyword{hplot}