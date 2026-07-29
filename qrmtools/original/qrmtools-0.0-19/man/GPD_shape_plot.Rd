\name{GPD_shape_plot}
\alias{GPD_shape_plot}
\title{Fitted GPD Shape as a Function of the Threshold}
\description{
  Fit GPDs to various thresholds and plot the fitted GPD shape as
  a function of the threshold.
}
\usage{
GPD_shape_plot(x, thresholds = seq(quantile(x, 0.5), quantile(x, 0.99),
                                   length.out = 65),
               estimate.cov = TRUE, conf.level = 0.95,
               CI.col = adjustcolor(1, alpha.f = 0.2),
               lines.args = list(), xlim = NULL, ylim = NULL,
               xlab = "Threshold", ylab = NULL,
               xlab2 = "Excesses", plot = TRUE, ...)
}
\arguments{
  \item{x}{\code{\link{vector}} of \code{\link{numeric}} data.}
  \item{thresholds}{\code{\link{numeric}} vector of thresholds for which
    to fit a GPD to the excesses.}
  \item{estimate.cov}{\code{\link{logical}} indicating whether
    confidence intervals are to be computed.}
  \item{conf.level}{confidence level of the confidence intervals if
    \code{estimate.cov}.}
  \item{CI.col}{color of the pointwise asymptotic confidence intervals
    (CIs); if \code{NA}, no CIs are shown.}
  \item{lines.args}{\code{\link{list}} of arguments passed to
    the underlying \code{\link{lines}()} for drawing the shape
    parameter as a function of the threshold.}
  \item{xlim, ylim, xlab, ylab}{see \code{\link{plot}()}.}
  \item{xlab2}{label of the secondary x-axis.}
  \item{plot}{\code{\link{logical}} indicating whether a plot is produced.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{plot}()}.}
}
\value{
  Invisibly returns a \code{\link{list}} containing the thresholds
  considered, the corresponding excesses and the fitted GPD
  objects as returned by the underlying \code{\link{fit_GPD_MLE}()}.
}
\details{
  Such plots can be used in the peaks-over-threshold method for
  determining the optimal threshold (as the smallest after which the plot
  is (roughly) stable).
}
\author{Marius Hofert}
\examples{
set.seed(271)
X <- rt(1000, df = 3.5)
GPD_shape_plot(X)
}
\keyword{distribution}
\keyword{hplot}