\name{Acf}
\alias{Acf}
\alias{plot.Acf}
\title{
Autocorrelation Function
}
\description{
  Plot the ACF without the traditional noninformation unit spike at lag
  0.
}
\usage{
Acf(x, lag.max = NULL, type = c("correlation", "covariance", "partial"),
    plot = TRUE, na.action = na.fail, demean = TRUE, \dots)

\method{plot}{Acf}(x, ci = 0.95, type = "h", xlab = "Lag", ylab = NULL, ylim = NULL, 
     main = NULL, ci.col = "blue", ci.type = c("white", "ma"), 
     max.mfrow = 6, ask = Npgs > 1 && dev.interactive(), 
     mar = if (nser > 2) c(3, 2, 2, 0.8) else par("mar"), 
     oma = if (nser > 2) c(1, 1.2, 1, 1) else par("oma"), 
     mgp = if (nser > 2) c(1.5, 0.6, 0) else par("mgp"), 
     xpd = par("xpd"), cex.main = if (nser > 2) 1 else par("cex.main"), 
     verbose = getOption("verbose"), acfLag0 = FALSE, \dots)
}
\arguments{
  \item{x}{
    for 'acf':  a numeric vector or time series.

    for 'plot.acf':  an object of class 'acf'.
  }
  \item{lag.max}{
    maximum lag at which to calculate the acf.
  }
  \item{ci}{
    coverage probability for confidence interval for 'plot.acf'.
  }
  \item{type}{
    the type of 'acf' or 'plot'
  }
  \item{plot}{
    logical. If 'TRUE' the 'acf' function will call 'plot.acf'.
  }
  \item{na.action}{
    function to be called by 'acf' to handle missing values.
  }
  \item{demean}{
    logical:  Should the x be replaced by \code{x - mean(x)} before computing
    the sums of squares and lagged cross products to produce the 'acf'?
  }
  \item{xlab,ylab,ylim,main,ci.col,ci.type,max.mfrow,ask,mar,oma,mgp,xpd,cex.main,verbose}{
    see the help page of \code{acf}: \code{help('acf', package = 'stats')}.
  }
  \item{acfLag0}{
    logical:  TRUE to plot the traditional noninformation unit spike at
    lag 0.  FALSE to omit that spike, consistent with the style in Tsay
    (2005).
  }
  \item{...}{further arguments passed to 'plot.acf'.}
}
\details{

  These functions are provided to make it easy to plot an
  autocorrelation function without the noninformative unit spike at lag 0. 
  This is done by calling plot(x, acfLag0 = FALSE, ...).  Apart from the
  'acfLag0' argument, the rest of the arguments are identical to those
  for 'acf' and 'plot.acf'.

}
\value{
  for \code{acf}, an object of class 'Acf', which inherits
  from class 'acf', as described with \code{help('acf', package='stats')}.

  for \code{plot.Acf}, NULL
}
\author{
  Spencer Graves for the FinTS modification of 'plot.acf'.
}
\references{
  Ruey Tsay (2005) Analysis of Financial Time Series, 2nd ed. (Wiley)
}
\seealso{
  \code{\link[stats]{acf}}
  \code{\link[stats]{plot.acf}}
  \code{\link{Box.test}}
  \code{\link{AutocorTest}}
}
\examples{
data(m.ibm2697)
Acf(m.ibm2697)
Acf(m.ibm2697, lag.max=100)
Acf(m.ibm2697, lag.max=100, main='Monthly IBM returns, 1926-1997')
}
\keyword{ts}
