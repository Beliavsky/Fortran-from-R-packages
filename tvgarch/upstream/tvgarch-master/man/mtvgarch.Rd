\name{mtvgarch}
\alias{mtvgarch}

\title{Estimate a multivariate TV-GARCH-X model}
 
\description{ 
Equation by equation estimation of a multivariate multiplicative TV-GARCH-X 
model with dnamic conditional correlations. For each variance equation, the 
long-term or unconditional component (TV) and the short-term or conditional 
variance component (GARCH-X) are estimated separately using maximization by 
parts, where the iterative algorithm proceeds until convergence. Conditional on 
the variance estimates, the dynamic conditional correlations are estimated by 
maximum likelihood.
}

\usage{
mtvgarch(y, order.g = c(1, 1), order.h = NULL, order.x = NULL, 
initial.values = list(), xtv = NULL, xreg = NULL, opt = 2, upper.speed = NULL,
tvgarch = FALSE, dcc = FALSE, turbo = TRUE, trace = FALSE)
}

\arguments{
\item{y}{numeric matrix, time series or \code{\link[zoo]{zoo}} object.}
\item{order.g}{integer matrix with each row indicating the order.g for each 
series; number of locations in each transition function of the TV components.}
\item{order.h}{integer matrix with each row indicating the order.h for each 
series; the first column controls the GARCH order, the second the ARCH order and 
the third the asymmetry order of the GARCH-X components. If \code{NULL}, the 
default, all series are assumed to follow a GARCH(1,1,0).}
\item{order.x}{\code{NULL} or binary matrix indicating which xreg variables 
should be included as covariates in the GARCH-X components. If provided and 
\code{xreg} is \code{NULL}, then the selected volatility spillovers are included 
as covariates.}
\item{initial.values}{a list containing the initial parameter values passed on 
to the optimisation routines (\link{constrOptim} for the TV component and 
\link{nlminb} for the GARCH-X component). If list(), the default, then the 
values are chosen automatically. TV component: \code{intercept.g} - \code{NULL} 
or numeric vector, \code{size} - \code{NULL} or numeric matrix containing the 
size initial coefficients, \code{speed} - \code{NULL} or numeric matrix 
containing the speed initial coefficients, \code{location} - \code{NULL} or 
numeric matrix containing the location initial coefficients. GARCH-X component: 
\code{intercept.h} - numeric vector, \code{arch} - \code{NULL} or numeric matrix 
containing the ARCH initial coefficients, \code{garch} - \code{NULL} or numeric 
matrix containing the GARCH-type initial coefficients, \code{asym} - \code{NULL} 
or numeric matrix containing the asymmetry-type initial coefficients, 
\code{par.xreg} - \code{NULL} or numeric matrix containing the X-type initial 
coefficients, and \code{R} - initial correlation coefficients.}
\item{xtv}{\code{NULL} or numeric vector, time series or zoo object to include 
as the transition variable in the TV component. If \code{NULL}, calendar time, 
scaled between 0 and 1, is used as the transition variable.}
\item{xreg}{numeric vector, time series or zoo object to include as covariates 
in the GARCH-X component.}
\item{opt}{integer indicating whether the speed parameter in the TV component 
should be scaled. If 0, no scaling; if 1, speed/sd(xtv); if 2, exp(speed).}
\item{upper.speed}{\code{NULL} or numeric \code{scalar} that sets 
the upper bound for \code{speed} in each transition function. If \code{NULL}, 
the default, the upper bound is 10000 for all transition functions. 
If numeric \code{scalar}, \code{upper.speed} is used for all transition 
functions.}
\item{tvgarch}{\code{logical}. If \code{TRUE}, the full parameter set
is estimated in one final step as well as the standard errors. If \code{FALSE} 
(default), estimates from last iteration are reported instead.}
\item{dcc}{logical. If \code{TRUE}, dynamic conditional correlations are 
estimated. If \code{FALSE} (default), then the conditional correlations are 
constant.}
\item{turbo}{logical. If \code{FALSE} (default), then the coefficient 
variance-covariance is computed during estimation, and the fitted values and 
residuals are attached to the returned object. If \code{TRUE}, then these 
operations are skipped, and hence estimation is faster. Note, however, that if 
turbo is set to \code{TRUE}, then the coefficient-covariance, fitted values and 
residuals can still be extracted subsequent to estimation with vcov.mtvgarch(), 
fitted.mtvgarch() and residuals.mtvgarch(), respectively.}
\item{trace}{logical. If \code{TRUE} all output is printed.}
}

\value{
An object of class 'mtvgarch'.
}

\references{
Cristina Amado and Timo Terasvirta (2013) Modelling volatility by variance 
decomposition, Journal of Econometrics 175, 142-153. 

Susana Campos-Martins and Genaro Sucarrat (2024) Modeling Nonstationary 
Financial Volatility with the R Package tvgarch, Journal of Statistical 
Software 108, 1-38.

Christian Francq and Jean-Michel Zakoian (2016) Estimating multivariate 
volatility models equation by equation, J. R. Stat. Soc. Ser. B Stat. Methodol 
78, 613-635.

Robert F. Engle (2002) Dynamic conditional correlation: A simple class of 
multivariate generalized autoregressive conditional heteroskedasticity models, 
Journal of Business and Economic Statistics 20, 339-350.
}

\author{Susana Campos-Martins}

\seealso{
\code{\link{tvgarch}},
\code{\link[garchx]{garchx}},
\code{\link{nlminb}},
\code{\link{constrOptim}}
}

\examples{
set.seed(12345)

## Simulate from a bivariate CCC-TV(1)-GARCH(1,1) model (default):
mySim <- mtvgarchSim(n = 1000)

## Estimate a CCC-TV(1)-GARCH(1,1) model (default):
myEst <- mtvgarch(y = mySim)

## Print estimation results:
print(myEst)

## Extract coefficients:
coef(myEst)

## Plot conditional volatilities:
plot(myEst)

## Generate predictions:
predict(myEst)
}

\keyword{Time Varying Parameter Models}
\keyword{Nonlinear Time Series}
\keyword{Financial Econometrics}
\keyword{Spillovers}

