\name{tvgarch}
\alias{tvgarch}

\title{Estimate a TV-GARCH-X model}
 
\description{ 
Quasi Maximum Likelihood (ML) estimation of a univariate multiplicative 
TV(s)-GARCH(p,q,r)-X model, where s indicates the number and the shape of the 
transition functions, r is the asymmetry order, p is the ARCH order, q is the 
GARCH order, and 'X' indicates that covariates can be included. Any transition 
variable, deterministic or stochastic, can be used to drive the transitions 
between the variance states. The TV long-term component introduces 
non-stationarity in the variance process, where the GARCH-X short-term component 
describes conditional heteroscedasticity. Maximization by parts leads to 
consistent and asymptotically normal estimates.
}

\usage{
tvgarch(y, order.g = 1, order.h = c(1,1,0), xtv = NULL, xreg = NULL, 
        initial.values = list(), opt = 2, upper.speed = NULL, tvgarch = FALSE,
        turbo = FALSE, trace = FALSE)
}

\arguments{
\item{y}{numeric \code{vector}, time series or \code{\link[zoo]{zoo}} object.}
\item{order.g}{\code{integer} \code{vector} of length s indicating the number of 
locations in each transition function of the TV component. Indicates whether a 
stationary GARCH or a nonstationary GARCH, i.e., TV-GARCH, shall be estimated.}
\item{order.h}{\code{integer} \code{vector} of the form c(p,q,r). The first 
entry controls the GARCH order, the second the ARCH order and the third the 
asymmetry order of the GARCH-X component.}
\item{initial.values}{a list containing the initial parameter values passed on 
to the optimisation routines (\link{constrOptim} for the TV component and 
\link{nlminb} for the GARCH-X component). If list(), the default, then the 
values are chosen automatically. TV component: \code{intercept.g} - \code{NULL} 
or numeric, \code{size} - \code{NULL} or numeric \code{vector} containing the 
size initial coefficients, \code{speed} - \code{NULL} or numeric \code{vector} 
containing the speed initial coefficients, \code{location} - \code{NULL} or 
numeric \code{vector} containing the location initial coefficients. GARCH-X or 
GARCH-X component of TV-GARCH-X: \code{intercept.h} - numeric, \code{arch} - 
\code{NULL} or numeric \code{vector} containing the ARCH initial coefficients, 
\code{garch} - \code{NULL} or numeric \code{vector} containing the GARCH-type 
initial coefficients, \code{asym} - \code{NULL} or numeric \code{vector} 
containing the assymetry-type initial coefficients, and \code{par.xreg} - 
\code{NULL} or numeric \code{vector} containing the X-type initial 
coefficients.}
\item{xtv}{\code{NULL} or numeric \code{vector}, time series or zoo object to 
include as the transition variable in the TV component. If \code{NULL}, calendar 
time, scaled between 0 and 1, is used as the transition variable. Not relevant 
for stationary GARCH.}
\item{xreg}{numeric \code{vector}, time series or zoo object to include as 
covariates in the GARCH-X component.}
\item{opt}{\code{integer} indicating whether the speed parameter in the TV 
component should be scaled. If 0, no scaling; if 1, \code{speed}/sd(\code{xtv}); 
if 2, exp(\code{speed}). Only relevant for TV-GARCH models.}
\item{upper.speed}{\code{NULL} or numeric \code{scalar} that sets 
the upper bound for \code{speed} in each transition function. If \code{NULL}, 
the default, the upper bound is 10000 for all transition functions. 
If numeric \code{scalar}, \code{upper.speed} is used for all transition 
functions.}
\item{tvgarch}{\code{logical}. If \code{TRUE}, the full parameter set
is estimated in one final step as well as the standard errors. If \code{FALSE} 
(default), estimates from last iteration are reported instead.}
\item{turbo}{\code{logical}. If \code{FALSE} (default), then the coefficient 
variance-covariance is computed during estimation, and the fitted values and 
residuals are attached to the returned object. If \code{TRUE}, then these 
operations are skipped, and hence estimation is faster. Note, however, that if 
turbo is set to \code{TRUE}, then the coefficient-covariance, fitted values and 
residuals can still be extracted subsequent to estimation with vcov.tvgarch(), 
fitted.tvgarch() and residuals.tvgarch(), respectively.}
\item{trace}{\code{logical}. If \code{TRUE} all output is printed when 
estimating a TV-GARCH.}
}

\value{
An object of class 'tvgarch'.
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
\code{\link[garchx]{garchx}},
\code{\link{tvgarchSim}},
\code{\link{nlminb}},
\code{\link{constrOptim}}
}

\examples{
set.seed(123)

## Simulate from a TV(1)-GARCH(1,1) model (default):
ySim <- tvgarchSim(n = 1500)

## Estimate a TV(1)-GARCH(1,1) model:
yEst <- tvgarch(y = ySim)

## Print estimation results:
print(yEst)

## Extract coefficients:
coef(yEst)

## Plot conditional volatilities:
plot(yEst)

## Extract log-likelihood:
logLik(yEst)

## Extract and store standardised residuals:
etaEst <- residuals(yEst)

## Generate predictions:
predict(yEst)
}

\keyword{Time Varying Parameter Models}
\keyword{Nonlinear Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
