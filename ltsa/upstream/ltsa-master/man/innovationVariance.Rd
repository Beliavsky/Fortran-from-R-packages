\name{innovationVariance}
\alias{innovationVariance}

\title{
Nonparametric estimate of the innovation variance
}
\description{
The innovation variance is estimated using a high order AR approximation determined by the AIC or
by using Kolmogoroff's formula with a smoothed periodogram. Default is AR.
}
\usage{
innovationVariance(z, method = c("AR", "Kolmogoroff"), ...)
}

\arguments{
  \item{z}{
time series
}
  \item{method}{
Default "AR". Set to "Kolmogoroff" for non-parametric periodogram estimate.
}
  \item{\dots}{
optional arguments that are passed to spec.pgram()
}
}
\value{ the innovation variance
}
\author{
A. I. McLeod
}

\seealso{
\code{\link{exactLoglikelihood}},
\code{\link{PredictionVariance}},
}

\examples{
z<-sunspot.year
#fitting high-order AR
innovationVariance(z)
#using periodogram
innovationVariance(z, method="Kolmogoroff")
#using smoothed periodogram
innovationVariance(z, method="Kolmogoroff", span=c(3, 3))
#the plot argument for spec.pgram() works too
innovationVariance(z, method="Kolmogoroff", span=c(3, 3), plot=TRUE)
}
\keyword{ ts }
