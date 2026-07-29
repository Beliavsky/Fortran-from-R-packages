\name{pp_qq_plot}
\alias{pp_plot}
\alias{qq_plot}
\title{P-P and Q-Q Plots}
\description{
  Probability-probability plots and quantile-quantile plots.
}
\usage{
pp_plot(x, FUN, pch = 20, xlab = "Theoretical probabilities",
        ylab = "Sample probabilities", ...)
qq_plot(x, FUN = qnorm, method = c("theoretical", "empirical"),
        pch = 20, do.qqline = TRUE, qqline.args = NULL,
        xlab = "Theoretical quantiles", ylab = "Sample quantiles",
        ...)
}
\arguments{
  \item{x}{data \code{\link{vector}}.}
  \item{FUN}{\code{\link{function}}. For
    \describe{
      \item{\code{pp_plot()}:}{The distribution function (vectorized).}
      \item{\code{qq_plot()}:}{The quantile function (vectorized).}
  }}
  \item{pch}{plot symbol.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{do.qqline}{\code{\link{logical}} indicating whether a Q-Q line
    is plotted.}
  \item{method}{method used to construct the Q-Q line. If
    \code{"theoretical"}, the theoretically true line with intercept 0
    and slope 1 is displayed; if \code{"empirical"}, the intercept
    and slope are determined with \code{\link{qqline}()}. The former
    helps deciding whether \code{x} comes from the distribution
    specified by \code{FUN} exactly, the latter whether \code{x}
    comes from a location-scale transformed distribution specified by
    \code{FUN}.}
  \item{qqline.args}{\code{\link{list}} containing additional arguments
    passed to the underlying \code{\link{abline}()} functions. Defaults to
    \code{list(a = 0, b = 1)} if \code{method = "theoretical"} and \code{list()}
    if \code{method = "empirical"}.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{plot}()}.}
}
\value{
  \code{\link{invisible}()}.
}
\details{
  Note that Q-Q plots are more widely used than P-P plots
  (as they highlight deviations in the tails more clearly).
}
\author{Marius Hofert}
\examples{
## Generate data
n <- 1000
mu <- 1
sig <- 3
nu <- 3.5
set.seed(271) # set seed
x <- mu + sig * sqrt((nu-2)/nu) * rt(n, df = nu) # sample from t_nu(mu, sig^2)

## P-P plot
pF <- function(q) pt((q - mu) / (sig * sqrt((nu-2)/nu)), df = nu)
pp_plot(x, FUN = pF)

## Q-Q plot
qF <- function(p) mu + sig * sqrt((nu-2)/nu) * qt(p, df = nu)
qq_plot(x, FUN = qF)

## A comparison with R's qqplot() and qqline()
qqplot(qF(ppoints(length(x))), x) # the same (except labels)
qqline(x, distribution = qF) # slightly different (since *estimated*)

## Difference of the two methods
set.seed(271)
z <- rnorm(1000)
## Standardized data
qq_plot(z, FUN = qnorm) # fine
qq_plot(z, FUN = qnorm, method = "empirical") # fine
## Location-scale transformed data
mu <- 3
sig <- 2
z. <- mu+sig*z
qq_plot(z., FUN = qnorm) # not fine (z. comes from N(mu, sig^2), not N(0,1))
qq_plot(z., FUN = qnorm, method = "empirical") # fine (as intercept and slope are estimated)
}
\keyword{hplot}