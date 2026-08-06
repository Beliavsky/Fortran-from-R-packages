\name{coef.tvgarchTest}

\alias{coef.tvgarchTest}
\alias{fitted.tvgarchTest}
\alias{logLik.tvgarchTest}
\alias{nobs.tvgarchTest}
\alias{plot.tvgarchTest}
\alias{predict.tvgarchTest}
\alias{print.tvgarchTest}
\alias{quantile.tvgarchTest}
\alias{residuals.tvgarchTest}
\alias{summary.tvgarchTest}
\alias{toLatex.tvgarchTest}
\alias{vcov.tvgarchTest}

\title{Extraction functions for univarate 'tvgarchTest' objects}

\description{Extraction functions for objects of class 'tvgarchTest'. Results 
from the estimation of the model under the null hypothesis, i.e., a GARCH(1,1) 
model, can be extracted similar to an object of class 'tvgarch' with the 
exception of functions print.tvgarchTest() and summary.tvgarchTest().}

\usage{
  \method{coef}{tvgarchTest}(object, ...)
  \method{fitted}{tvgarchTest}(object, 
  as.zoo = TRUE, ...)
  \method{logLik}{tvgarchTest}(object, ...)
  \method{nobs}{tvgarchTest}(object, ...)
  \method{plot}{tvgarchTest}(x, ...)                           
  \method{predict}{tvgarchTest}(object, n.ahead = 10, newxreg = NULL, 
                                 newindex = NULL, n.sim = 5000,
                                 as.zoo = TRUE, verbose = FALSE, ...)
  \method{print}{tvgarchTest}(x, ...)
  \method{quantile}{tvgarchTest}(x, probs = 0.025, names = TRUE, type = 7, 
                                  as.zoo = TRUE, ...)
  \method{residuals}{tvgarchTest}(object, as.zoo = TRUE, ...)
  \method{summary}{tvgarchTest}(object, ...)
  \method{toLatex}{tvgarchTest}(object, digits = 4, ...)
  \method{vcov}{tvgarchTest}(object, ...)
}

\arguments{
  \item{object}{an object of class 'tvgarchTest'.}
  \item{x}{an object of class 'tvgarchTest'.}
  \item{as.zoo}{logical. If \code{TRUE}, then the returned result is of 
  class \code{\link[zoo]{zoo}}.}
  \item{n.ahead}{integer that determines how many steps ahead predictions 
  should be generated.}
  \item{newxreg}{vector or matrix with the out-of-sample regressor 
  values.}
  \item{newindex}{a \code{\link[zoo]{zoo}}-index for the out-of-sample predictions. If 
  \code{NULL} (default), then \code{1:n.ahead} is used.}
  \item{n.sim}{integer, the number of simulations.}
  \item{verbose}{logical. If \code{TRUE}, then the simulations - in 
  addition to the predictions - are returned.}
  \item{probs}{vector of probabilities.}
  \item{names}{logical, whether to return names or not.}
  \item{type}{integer that determines the algorithm used to compute the 
  quantile, see \code{\link{quantile}}.}
  \item{digits}{integer, the number of digits in the printed LaTeX code.}
  \item{\dots}{additional arguments.}
}

\value{
  \item{coef:}{parameter estimates.}
  \item{fitted:}{fitted conditional variance.}
  \item{logLik:}{optimised log-likelihood (normal density) value.}
  \item{nobs:}{the number of observations used in the estimation.}
  \item{plot:}{plot of the fitted conditional volatility.}
  \item{predict:}{variance predictions.}
  \item{quantile:}{fitted quantiles, i.e. the conditional standard deviation 
  times the empirical quantile of the standardised innovations.}
  \item{residuals:}{volatility standardised residuals.}
  \item{summary:}{summary of test result.}
  \item{vcov:}{coefficient variance-covariance matrix.}
}

\references{
Cristina Amado and Timo Terasvirta (2013) Modelling volatility by variance 
decomposition, Journal of Econometrics 175, 142-153. 

Cristina Amado and Timo Terasvirta (2014) Modelling changes in the 
unconditional variance of long stock return series, Journal of Empirical Finance 
25, 15-35. 

Susana Campos-Martins and Genaro Sucarrat (2024) Modeling Nonstationary 
Financial Volatility with the R Package tvgarch, Journal of Statistical 
Software 108, 1-38.
}

\author{Susana Campos-Martins}

\seealso{
\code{\link{tvgarchTest}}, 
\code{\link{tvgarchSim}}, 
\code{\link{tvgarch}}, 
\code{\link[zoo]{zoo}}
}

\examples{
set.seed(123)

## Simulate from a TV(1)-GARCH(1,1) model (default):
ySim <- tvgarchSim(n = 1500)

## Test a GARCH(1,1) model against a TV(1)-GARCH(1,1) model:
yTest <- tvgarchTest(y = ySim)

## Print test and estimation results:
print(yTest)

## Estimated number of locations
summary(yTest)

## Extract and plot estimation results for GARCH(1,1) used in the test:
sigma2Test <- fitted(yTest)
plot(yTest)

## Estimate a TV(s)-GARCH(1,1) model:
s <- summary(yTest)
yEst <- tvgarch(y = ySim, order.g  = s)
}

\keyword{Time Varying Parameter Models}
\keyword{Nonlinear Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
