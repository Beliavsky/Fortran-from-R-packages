\name{exactLoglikelihood}
\alias{exactLoglikelihood}

\title{Exact log-likelihood and MLE for variance}

\description{
Provides an exact log-likelihood that is exactly equal to the value of the
probability density function with the random variables replaced by data and the
parameters replaced by their estimated value. The corresponding estimate of the
variance term is return.
}
\usage{
exactLoglikelihood(r, z, innovationVarianceQ = TRUE)
}

\arguments{
  \item{r}{the portion of autocovariance function which when multiplied by
  the variance term equals the full autocovariance function.
}
  \item{z}{the time series assumed to have mean zero
}
  \item{innovationVarianceQ}{
When TRUE, the variance term is the innovation variance and when FALSE it is the
variance of the time series. For ARFIMA models, set to TRUE. But FGN requires
setting innovationVarianceQ to FALSE since only the innovation variance is not
known and so the likelihood has a slightly different form.
}
}
\details{
This function uses the trench algorithm that is implememented in C.
This function is provided to include all multiplicative constants.
For many purposes, such as MLE, we only need to likelihood function up to
a multiplicative constant. But for information criteria, we may need the constant
terms so we can compare our results with other types of models or with other software
such as arima(). The arima() function also computes the exact log-likelihood
and uses it in the computation of the AIC and BIC.
}
\value{

\item{LL }{exact log-likelihood}
\item{sigmaSq }{MLE for the variance term. If innovationVarianceQ is TRUE, is the an
estimate of the residual variance otherwise it is an estimate of the variance of the time series.}
}

\author{
A. I. McLeod, aimcleod@uwo.ca
}

\seealso{
\code{\link{TrenchLoglikelihood}},
\code{\link{DLLoglikelihood}}
}
\examples{
set.seed(7773311)
n <- 200
z <- arima.sim(model=list(ar=0.9, ma=-0.6), n=n, n.start=10^4)
out <- arima(z, order=c(1,0,1), include.mean=FALSE)
out
#note
#sigma^2 estimated as 0.9558:  log likelihood = -279.66,  aic = 565.31
r <- tacvfARMA(phi=coef(out)[1], theta=-coef(out)[2], maxLag=n-1)
exactLoglikelihood(r, z, innovationVarianceQ = TRUE)
#agrees!
}
\keyword{ ts }
