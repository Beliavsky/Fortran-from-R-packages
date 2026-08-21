\name{ARMAspec}
\alias{ARMAspec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Theoretical spectral density function of a stationary ARMA model}
\description{Computes and plots the theoretical spectral density function of a stationary ARMA model}
\usage{
ARMAspec(model, freq = seq(0, 0.5, 0.001), plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{an arma model}
  \item{freq}{  vector of frequency over which the spectral density is computed}
  \item{plot}{if true, plot the spectral density function; default is true}
  \item{\dots}{other parameters to be passed to the plot function}
}
\value{a list:
  \item{spec}{spectral density values}
  \item{freq}{same as freq in the input}
  \item{model}{the arma model}
}
\author{Kung-Sik Chan}
\seealso{\code{\link{spec}}}
\examples{
\donttest{
theta=.9 # Reset theta for other MA(1) plots
ARMAspec(model=list(ma=-theta)) 
}
}
\keyword{methods}
