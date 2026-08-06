\name{coef.tvgarch}

\alias{coef.tvgarch}
\alias{fitted.tvgarch}
\alias{logLik.tvgarch}
\alias{nobs.tvgarch}
\alias{plot.tvgarch}
\alias{predict.tvgarch}
\alias{print.tvgarch}
\alias{quantile.tvgarch}
\alias{residuals.tvgarch}
\alias{summary.tvgarch}
\alias{toLatex.tvgarch}
\alias{vcov.tvgarch}

\title{Extraction functions (S3 methods) for univarate 'tvgarch' objects}

\description{Extraction functions (S3 methods) for objects of class 'tvgarch'.}

\usage{
  \method{coef}{tvgarch}(object, spec = c("tvgarch", "garch", "tv"), ...)
  \method{fitted}{tvgarch}(object, spec = c("tvgarch", "garch", "tv"), 
                            as.zoo = TRUE, ...)
  \method{logLik}{tvgarch}(object, ...)
  \method{nobs}{tvgarch}(object, ...)
  \method{plot}{tvgarch}(x, spec = c("tvgarch", "garch", "tv"), ...)                            
  \method{predict}{tvgarch}(object, n.ahead = 10, newxtv = NULL, 
                             newxreg = NULL, newindex = NULL, n.sim = 5000, 
                             as.zoo = TRUE, verbose = FALSE, ...)
  \method{print}{tvgarch}(x, ...)
  \method{quantile}{tvgarch}(x, probs = 0.025, names = TRUE, type = 7, 
                              as.zoo = TRUE, ...)
  \method{residuals}{tvgarch}(object, as.zoo = TRUE, ...)
  \method{summary}{tvgarch}(object, ...)
  \method{toLatex}{tvgarch}(object, digits = 4, ...)
  \method{vcov}{tvgarch}(object, spec = c("tvgarch", "garch", "tv"), ...)
}

\arguments{
  \item{object}{an object of class 'tvgarch'.}
  \item{spec}{specifies whether the function should extract specific results. If 
  "tv", extracts results for the TV component and if "garch" extracts results 
  for the GARCH-X component of TV-GARCH-X model. If "tvgarch", extracts 
  results for TV-GARCH-X model. Only relevant for TV-GARCH-X models. 
  Otherwise, extracts results for GARCH-X models.}
  \item{x}{an object of class 'tvgarch'.}
  \item{as.zoo}{logical. If \code{TRUE}, then the returned result is of 
  class \code{\link[zoo]{zoo}}.}
  \item{n.ahead}{integer that determines how many steps ahead predictions 
  should be generated.}
  \item{newxtv}{\code{NULL} or vector with the out-of-sample transition 
  variable. If \code{NULL}, out-of-sample component g equals \code{intercept.g}. 
  Only relevant for TV-GARCH-X models.}
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
  \item{print:}{print of the estimation results.}
  \item{quantile:}{fitted quantiles, i.e. the conditional standard deviation 
  times the empirical quantile of the standardised innovations.}
  \item{residuals:}{volatility standardised residuals.}
  \item{summary:}{summary of estimation results.}
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
\code{\link{tvgarch}}, 
\code{\link{tvgarchSim}}, 
\code{\link[zoo]{zoo}}
}

\examples{
set.seed(123)

## Simulate from a TV(1)-GARCH(1,1) model (default):
ySim <- tvgarchSim(n = 1500)

## Estimate a TV(1)-GARCH(1,1) model:
yEst <- tvgarch(y = ySim)

## Print estimation results:
print(yEst)

## Extract and store conditional variances:
sigma2Est <- fitted(yEst)

## Plot:
plot(yEst)

## Generate predictions:
predict(yEst)
}

\keyword{Time Varying Parameter Models}
\keyword{Nonlinear Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
