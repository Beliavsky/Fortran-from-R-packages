\name{coef.mtvgarch}

\alias{coef.mtvgarch}
\alias{fitted.mtvgarch}
\alias{logLik.mtvgarch}
\alias{nobs.mtvgarch}
\alias{plot.mtvgarch}
\alias{predict.mtvgarch}
\alias{print.mtvgarch}
\alias{quantile.mtvgarch}
\alias{residuals.mtvgarch}
\alias{summary.mtvgarch}
\alias{toLatex.mtvgarch}
\alias{vcov.mtvgarch}

\title{Extraction functions for multivariate 'mtvgarch' objects}

\description{Extraction functions for objects of class 'mtvgarch'.}

\usage{
  \method{coef}{mtvgarch}(object, spec = c("tvgarch", "garch", "tv", "cc"), ...)
  \method{fitted}{mtvgarch}(object, spec = c("tvgarch", "garch", "tv", "cc"), 
                             as.zoo = TRUE, ...)
  \method{logLik}{mtvgarch}(object, ...)
  \method{nobs}{mtvgarch}(object, ...)
  \method{plot}{mtvgarch}(x, spec = c("tvgarch", "garch", "tv"), ...)                           
  \method{predict}{mtvgarch}(object, n.ahead = 10, newxtv = NULL, 
                              newxreg = NULL, newindex = NULL, n.sim = 5000, 
                              as.zoo = TRUE, verbose = FALSE, ...)
  \method{print}{mtvgarch}(x, ...)
  \method{quantile}{mtvgarch}(x, probs = 0.025, type = 7, as.zoo = TRUE, ...)
  \method{residuals}{mtvgarch}(object, as.zoo = TRUE, ...)
  \method{summary}{mtvgarch}(object, ...)
  \method{toLatex}{mtvgarch}(object, digits = 4, ...)
  \method{vcov}{mtvgarch}(object, spec = c("tvgarch", "garch", "tv", "cc"), ...)
}

\arguments{
  \item{object}{an object of class 'mtvgarch'.}
  \item{spec}{specifies whether the function should extract specific results. If 
  "tv", extracts results for the TV component and if "garch" extracts results 
  for the GARCH-X component of TV-GARCH-X model. If "tvgarch", extracts 
  results for TV-GARCH-X model. Only relevant for TV-GARCH-X models. 
  Otherwise, extracts results for GARCH-X models. If "cc", extracts results 
  concerning the conditional correlations. "cc" is not valid in plot.mtvgarch().}
  \item{x}{an object of class 'mtvgarch'.}
  \item{as.zoo}{logical. If \code{TRUE}, then the returned result is of class 
  \code{\link[zoo]{zoo}}.}
  \item{n.ahead}{integer that determines how many steps ahead predictions 
  should be generated.}
  \item{newxtv}{\code{NULL} or vector with the out-of-sample transition 
  variable. If \code{NULL}, out-of-sample component g equals \code{intercept.g}. 
  Only relevant for TV-GARCH-X models.}
  \item{newxreg}{vector or matrix with the out-of-sample regressor values.}
  \item{newindex}{\code{\link[zoo]{zoo}}-index for the out-of-sample predictions. If 
  \code{NULL} (default), then \code{1:n.ahead} is used.}
  \item{n.sim}{integer, the number of simulations.}
  \item{verbose}{logical. If \code{TRUE}, then the simulations - in 
  addition to the predictions - are returned.}
  \item{probs}{vector of probabilities.}
  \item{type}{integer that determines the algorithm used to compute the 
  quantile, see \code{\link{quantile}}.}
  \item{digits}{integer, the number of digits in the printed LaTeX code.}
  \item{\dots}{additional arguments.}
}

\value{
  \item{coef:}{parameter estimates.}
  \item{fitted:}{fitted conditional variances and correlations.}
  \item{logLik:}{optimised log-likelihood (normal density) values.}
  \item{nobs:}{number of observations used in the estimation.}
  \item{plot:}{plots of the fitted conditional volatilities.}
  \item{predict:}{variance predictions. Column order differs when spillovers are 
  allowed.}
  \item{print:}{print of the estimation results.}
  \item{quantile:}{fitted quantiles, i.e. the conditional standard deviation 
  times the empirical quantile of the standardised innovations.}
  \item{residuals:}{volatility standardised residuals.}
  \item{summary:}{summary of estimation results.}
  \item{vcov:}{coefficient variance-covariance matrices.}
}

\references{
Cristina Amado and Timo Terasvirta (2013) Modelling volatility by variance 
decomposition, Journal of Econometrics 175, 142-153. 

Christian Francq and Jean-Michel Zakoian (2016) Estimating multivariate 
volatility models equation by equation, J. R. Stat. Soc. Ser. B Stat. Methodol 
78, 613-635.

Susana Campos-Martins and Genaro Sucarrat (2024) Modeling Nonstationary 
Financial Volatility with the R Package tvgarch, Journal of Statistical 
Software 108, 1-38.
}

\author{Susana Campos-Martins}

\seealso{
\code{\link{mtvgarch}}, 
\code{\link{mtvgarchSim}}, 
\code{\link{tvgarch}}, 
\code{\link[garchx]{garchx}}, 
\code{\link[zoo]{zoo}}
}

\examples{
set.seed(12345)

## Simulate from a bivariate CCC-TV(1)-GARCH(1,1) model (default):
mySim <- mtvgarchSim(n = 1500)

## Estimate a CCC-TV(1)-GARCH(1,1) model:
myEst <- mtvgarch(y = mySim)

## Print estimation results:
print(myEst)

## Extract and store conditional variances:
sigma2Est <- fitted(myEst)

## Plot:
plot(myEst)

## Generate predictions:
predict(myEst)
}

\keyword{Time Varying Parameter Models}
\keyword{Nonlinear Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
