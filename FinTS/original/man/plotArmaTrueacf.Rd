\name{plotArmaTrueacf}
\alias{plotArmaTrueacf}
\title{ Plot the theoretical ACF corresponding to an ARMA model  }
\description{
  Compute the roots and theoretical ACF corresponding to an ARMA model
}
\usage{
plotArmaTrueacf(object, lag.max = 20, pacf = FALSE, plot = TRUE,
          xlab = "lag", ylab = c("ACF", "PACF")[1+pacf],
          ylim = c(-1, 1) * max(ACF), type = "h",
          complex.eps = 1000 * .Machine[["double.neg.eps"]], ...)
}
\arguments{
  \item{object}{
    either a numeric vector or a list with components 'ar' and 'ma'.  If
    'object' is numeric, it is interpreted as a model with no 'ma' part.
  }
  \item{lag.max}{
    the maximum number of lags for which to calculate the ACF or PACF.
  }
  \item{pacf}{
    logical.  Should the partial autocorrelations be returned?
  }
  \item{plot}{
    logical.  Should the ACF (or PACF) be plotted?
  }
  \item{xlab, ylab, ylim, type}{
    arguments for 'plot'
  }
  \item{complex.eps}{
    a small positive number used to identify complex conjugates:  Let
    'roots' = the vector of p roots of the characteristic polynomial of
    the autoregressive part of 'object'.  This is used by
    'findConjugates':  x[i] and x[j] are considered conjugate if their
    relative difference exceeds complex.eps but the relative difference
    of their conjugates is less than complex.eps.

    We use 'solve' in the 'polynom' package, because it was substially
    more accurate for cases we tested in R 2.6.0 than 'polyroot'.
  }
  \item{...}{
    optional arguments passed to 'plot'.
  }
}
\details{
  1.  Compute and test stationarity.  An ARMA process is stationary if
  all the roots of its AR component lie inside the unit  circle (Box and
  Jenkins, 1970).  If the process is not stationary, a warning is
  issued, and no plot is produced.

  2.  Compute and plot the theoretical ACF.

  3.  Analyze periodicity of any complex roots
}
\value{
  a list with the following components

  \item{roots}{
    a complex vector of the roots sorted by modulus and sign of the
    immaginary part.
  }
  \item{acf, pacf}{
    a named numeric vector of the estimated ACF (or PACF of 'pacf =
    TRUE').
  }
  \item{periodicity}{
    a data.frame with one row for each complex conjugate pair of roots
    and columns 'damping' and 'period'.
  }
}
\source{
  \url{https://faculty.chicagobooth.edu/ruey-s-tsay/teaching}
}
\references{
  George E. P. Box and Gwilym M. Jenkins (1970) Time Series Analysis,
  Forecasting and Control (Holden-Day, sec. 3.4.1.  Stationarity and
  invertibility properties of Mixed autoregressive-moving average
  processes)

  Ruey Tsay (2005) Analysis of Financial Time Series, 2nd ed. (Wiley,
  ch. 2)
}
\seealso{
  \code{\link[polynom]{solve.polynomial}},
  \code{\link[stats]{ARMAacf}}
}
\examples{
# Tsay, Figure 2.3
op <- par(mfcol = c(1, 2))
plotArmaTrueacf(.8, lag.max = 8)
title("(a)")
plotArmaTrueacf(-.8, lag.max = 8)
title("b")
par(op)

# Tsay, Figure 2.4
op <- par(mfrow = c(2,2))
plotArmaTrueacf(c(1.2, -.35))
title("(a)")
plotArmaTrueacf(c(.6, -.4))
title("(b)")
plotArmaTrueacf(c(.2, .35))
title("(c)")
plotArmaTrueacf(c(-.2, .35))
title("(d)")
par(op)

# Tsay, Example 2.1
data(q.gnp4791)
(fit.ar3 <- ar(q.gnp4791, aic = FALSE, order = 3))
plotArmaTrueacf(fit.ar3)
}
\keyword{hplot}
