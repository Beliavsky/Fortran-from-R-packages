\name{tail_plot}
\alias{tail_plot}
\title{Plot of an Empirical Surival Function with Smith Estimator}
\description{
  Plot an empirical tail survival function, possibly overlaid
  with the Smith estimator.
}
\usage{
tail_plot(x, threshold, shape = NULL, scale = NULL,
          q = NULL, length.out = 129, lines.args = list(),
          log = "xy", xlim = NULL, ylim = NULL,
          xlab = "x", ylab = "Tail probability at x", ...)
}
\arguments{
  \item{x}{\code{\link{numeric}} vector of data.}
  \item{threshold}{\code{numeric(1)} giving the threshold \eqn{u} above
   which the tail (starts and) is to be plotted.}
  \item{shape}{\code{NULL} or the GPD shape parameter \eqn{\xi}{xi}
    (typically obtained via \code{\link{fit_GPD_MLE}()}).}
  \item{scale}{\code{NULL} or the GPD shape parameter \eqn{\beta}{beta}
    (typically obtained via \code{\link{fit_GPD_MLE}()}).}
  \item{q}{\code{NULL}, \code{numeric(1)} or \code{\link{numeric}} vector of
    evaluationn points of the  Smith estimator (semi-parametric GPD-based tail
    estimator in the POT method). If \code{NULL}, the evaluation points
    are determined internally as an equidistant sequence of length
    \code{length.out} between the smallest and largest exceedance
    (taken equidistant in log-scale if \code{log} contains \code{"x"}).
    If \code{numeric(1)}, then the behavior is similar to \code{NULL}
    with the exception that the plot is extended to the right of the
    largest exceedance if \code{q} is larger than the largest
    exceedance.}
  \item{length.out}{length of \code{q}.}
  \item{lines.args}{\code{\link{list}} of arguments passed to the
    underlying \code{\link{lines}()}.}
  \item{log}{\code{\link{character}} indicating whether
    logarithmic axes are to be used.}
  \item{xlim}{x-axis limits.}
  \item{ylim}{y-axis limits.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{plot}()}.}
}
\value{
  If both \code{shape} and \code{scale} are provided, \code{tail_plot()}
  overlays the empirical tail survival function estimator (evaluated at
  the exceedances) with the corresponding GPD. In this case,
  \code{tail_plot()} invisibly returns a list with two two-column
  matrices, one containing the x-values and y-values of the
  empirical survival distribution estimator and one containing the
  x-values and y-values of the Smith estimator. If \code{shape} or
  \code{scale} are \code{NULL}, \code{tail_plot()} invisibly returns
  a two-column matrix with the x-values and y-values of the empirical
  survival distribution estimator.
}
\author{Marius Hofert}
\examples{
## Generate losses to work with
set.seed(271)
X <- rt(1000, df = 3.5) # in MDA(H_{1/df}); see MFE (2015, Section 16.1.1)

## Threshold (see ?dGPDtail, for example)
u <- 1.5 # threshold

## Plots of empirical survival distribution functions (overlaid with Smith estimator)
tail_plot(X, threshold = u, log = "", type = "b") # => need log-scale
tail_plot(X, threshold = u, type = "s") # as a step function
fit <- fit_GPD_MLE(X[X > u] - u) # fit GPD to excesses (POT method)
tail_plot(X, threshold = u, # without log-scale
          shape = fit$par[["shape"]], scale = fit$par[["scale"]], log = "")
tail_plot(X, threshold = u, # highlights linearity
          shape = fit$par[["shape"]], scale = fit$par[["scale"]])
}
\keyword{hplot}