\name{BoxCox.ar}
\alias{BoxCox.ar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Determine the power transformation for serially correlated data}
\description{
Determine the appropriate 
 power transformation for time-series data. The objective is
to estimate the power transformation so that the transformed time series is
approximately a Gaussian AR process.} 
\usage{
BoxCox.ar(y, order, lambda = seq(-2, 2, 0.01), plotit = TRUE, 
method = c("mle", "yule-walker", "burg", "ols", "yw"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ univariate time series (must be positive)}
  \item{order}{ AR order for the data; if missing, the order is determined by
   AIC for the log-transformed data}
  \item{lambda}{ a vector of candidate power transformation values; if missing, it is
set to be from -2 to 2, with increment .01}
  \item{plotit}{ logical value, if true, plot the profile
log-likelihood for  
the power estimator}
  \item{method}{ method of AR estimation; default is "mle"}
  \item{\dots}{ other parameters to be passed to the ar function}
}
\value{A list that contains the following:
 \item{lambda}{candidate power transformation parameter values}
  \item{loglike}{profile log-likelihood}
  \item{mle}{maximum likelihood estimate of the power transformation
value}
  \item{ci}{95\% C.I. of the power transformation value} }
\author{Kung-Sik Chan}
\note{ The procedure is very computer intensive. Be patient for the outcome}
%\seealso{ \code{\link{boxcox}}}
\examples{
data(hare)
# hare.transf=BoxCox.ar(y=hare)
# hare.transf$ci
}
\keyword{methods}

