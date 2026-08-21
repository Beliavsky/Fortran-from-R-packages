\name{prewhiten}
\alias{prewhiten}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prewhiten a Bivariate Time Series, and Compute and Plot Their Sample Cross-Correlation Function}
\description{
The bivariate time series are prewhitened according to an AR model fitted to
the x-component of the bivariate series. Alternatively, if an ARIMA model 
 is provided, it will be used to prewhiten both series. The CCF of the
prewhitened bivariate
series is then computed and plotted.
}
\usage{
prewhiten(x, y, x.model = ar.res,ylab="CCF", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{first component series}
  \item{y}{second  component series}
  \item{x.model}{an ARIMA model; if provided, it is used to prewhiten both series. Otherwise, an AR model is fitted to the x-series and used to prewhiten both
series.The AR order is chosen by minimizing the AIC and the fit carried out by
the ar.ols function.}
  \item{ylab}{label of y-axis; default is "CCF"}
  \item{\dots}{ additional parameters to be passed to the ar.ols and the
ccf function.}
}
\value{A list containing the following components: 
  \item{ccf}{Output from the ccf function on the prewhitened data.}
  \item{ar}{The AR model fit to the x-series, or x.model if it is provided.}
}
\author{Kung-Sik Chan}
\examples{
data(milk)
data(electricity)
milk.electricity=ts.intersect(milk,log(electricity))
plot(milk.electricity,yax.flip=TRUE,main='') 
ccf(as.numeric(milk.electricity[,1]),as.numeric(milk.electricity[,2]),
main='milk & electricity',ylab='CCF')
me.dif=ts.intersect(diff(diff(milk,12)),diff(diff(log(electricity),12)))
prewhiten(as.numeric(me.dif[,1]),as.numeric(me.dif[,2]),
,ylab='CCF' )
}
\keyword{methods}
