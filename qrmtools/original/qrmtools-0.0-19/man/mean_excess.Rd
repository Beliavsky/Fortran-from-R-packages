\name{mean_excess}
\alias{mean_excess_np}
\alias{mean_excess_plot}
\alias{mean_excess_GPD}
\title{Mean Excess}
\description{
  Sample mean excess function, mean excess function of a GPD
  and sample mean excess plot.
}
\usage{
mean_excess_np(x, omit = 3)
mean_excess_plot(x, omit = 3,
                 xlab = "Threshold", ylab = "Mean excess over threshold", ...)
mean_excess_GPD(x, shape, scale)
}
\arguments{
  \item{x}{
    \describe{
      \item{\code{mean_excess_GPD()}:}{\code{\link{numeric}} vector of
	evaluation points of the mean excess function of the GPD.}
      \item{otherwise:}{\code{\link{numeric}} vector of data.}
    }
  }
  \item{omit}{number \eqn{\ge 1}{>= 1} of unique % differs from QRM::MEplot()!
    last observations to be omitted from the sorted data (as mean excess plot
    becomes unreliable for these observations as thresholds).}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{plot}()}.}
  \item{shape}{GPD shape parameter \eqn{\xi}{xi}.}
  \item{scale}{GPD scale parameter \eqn{\beta}{beta}.}
}
\value{
  \code{mean_excess_np()} returns a two-column matrix giving
  the sorted data without the \code{omit}-largest unique values
  (first column) and the corresponding values of the sample mean excess
  function (second column). It is mainly used in \code{mean_excess_plot()}.

  \code{mean_excess_plot()} returns \code{invisible()}.

  \code{mean_excess_GPD()} returns the mean excess function of a
  generalized Pareto distribution evaluated at \code{x}.
}
\details{
  Mean excess plots can be used in the peaks-over-threshold method for
  choosing a threshold. To this end, one chooses the smallest threshold
  above which the mean excess plot is roughly linear.
}
\author{Marius Hofert}
\examples{
## Generate losses to work with
set.seed(271)
X <- rt(1000, df = 3.5) # in MDA(H_{1/df}); see MFE (2015, Section 16.1.1)

## (Sample) mean excess plot and threshold choice
mean_excess_plot(X[X > 0]) # we only use positive values here to see 'more'
## => Any value in [0.8, 2] seems reasonable as threshold at first sight
##    but 0.8 to 1 turns out to be too small for the degrees of
##    freedom implied by the GPD estimator to be close to the true value 3.5.
## => We go with threshold 1.5 here.
u <- 1.5 # thresholds

## An alternative way
ME <- mean_excess_np(X[X > 0])
plot(ME, xlab = "Threshold", ylab = "Mean excess over threshold")

## Mean excess plot with mean excess function of the fitted GPD
fit <- fit_GPD_MLE(X[X > u] - u)
q <- seq(u, ME[nrow(ME),"x"], length.out = 129)
MEF.GPD <- mean_excess_GPD(q-u, shape = fit$par[["shape"]], scale = fit$par[["scale"]])
mean_excess_plot(X[X > 0]) # mean excess plot for positive losses...
lines(q, MEF.GPD, col = "royalblue", lwd = 1.4) # ... with mean excess function of the fitted GPD
}
\keyword{nonparametric}
\keyword{hplot}