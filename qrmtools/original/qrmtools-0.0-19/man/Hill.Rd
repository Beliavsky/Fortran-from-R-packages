\name{Hill}
\alias{Hill_estimator}
\alias{Hill_plot}
\title{Hill Estimator and Plot}
\description{
  Compute the Hill estimator and Hill plot.
}
\usage{
Hill_estimator(x, k = c(10, length(x)), conf.level = 0.95)
Hill_plot(x, k = c(10, length(x)), conf.level = 0.95, Hill.estimator = NULL,
          log = "x", xlim = NULL, ylim = NULL,
          xlab = "Order statistics", ylab = "Tail index",
          CI.col = adjustcolor(1, alpha.f = 0.2), lines.args = list(),
          xaxis2 = TRUE, xlab2 = "Empirical probability", ...)
}
\arguments{
  \item{x}{\code{\link{vector}} of \code{\link{numeric}} data.}
  \item{k}{\code{\link{vector}} of length 2, determining the smallest
    and largest number of order statistics of \code{x} to compute the
    Hill estimator for (the smallest needs to be >= 2). If \code{k} is
    of length 1, \code{k} is expanded by \code{length(x)}.}
  \item{conf.level}{confidence level of the confidence intervals.}
  \item{Hill.estimator}{object as returned by \code{Hill_estimator()}.}
  \item{log, xlim, ylim, xlab, ylab}{see \code{\link{plot}()}.}
  \item{CI.col}{color of the pointwise asymptotic confidence intervals
    (CIs); if \code{NA}, no CIs are shown.}
  \item{lines.args}{\code{\link{list}} of additional arguments
    for the underlying \code{\link{lines}()} call to draw the Hill
    estimator.}
  \item{xaxis2}{\code{\link{logical}} indicating whether a third
    axis is drawn that shows the empirical probabilities
    \code{1-(k-1)/length(x)} corresponding to \code{k}, so the value of
    the column \code{k.prob} as returned by \code{Hill_estimator()}.}
  \item{xlab2}{label of the secondary x-axis.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{plot}()}.}
}
\value{
  \describe{
    \item{\code{Hill_estimator()}:}{A five-column matrix containing the
      indices \code{k}, their corresponding empirical probabilities
      \code{k.prob}, the estimated tail indices \code{tail.index},
      and the lower and upper CI endpoints \code{CI.low} and \code{CI.up}.}
    \item{\code{Hill_plot()}:}{Hill plot by side-effect.}
  }
}
\details{
  See McNeil et al. (2015, Section 5.2.4, (5.23))
}
\author{Marius Hofert}
\references{
   McNeil, A. J., Frey, R. and Embrechts, P. (2015).
   \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
   Princeton University Press.
}
\examples{
set.seed(271)
X <- rt(1000, df = 3.5)
Y <- X[X > 0]
Hill_plot(Y)
Hill_plot(Y, log = "", CI.col = NA)
}
\keyword{nonparametric}
\keyword{hplot}