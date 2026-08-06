\name{tvgarchTest}
\alias{tvgarchTest}

\title{Test of a multiplicative time-varying GARCH model}

\description{  
Compute the non-robust and robust Lagrange-Multiplier (LM-)type test statistics 
for examining the null hypothesis of constant long-term variance, GARCH(1,1), 
against the alternative of a smoothly changing long-term component, 
TV-GARCH(1,1).
}

\usage{
tvgarchTest(y, xtv = NULL, alpha = 0.05)
}

\arguments{
\item{y}{numeric \code{vector}, time series or \code{\link[zoo]{zoo}} object.}
\item{xtv}{\code{NULL} or numeric vector, time series or zoo object to include 
as the transition variable in the TV component. If \code{NULL}, calendar time 
scaled between 0 and 1 is used as the transition variable.}
\item{alpha}{the significance level.}
}

\value{
An object of class 'tvgarchTest'.
}

\references{
Cristina Amado and Timo Terasvirta (2017) Specification and 
testing of multiplicative time varying GARCH models with applications, 
Econometric Reviews 36:4, 421-446.
}

\author{Susana Campos-Martins}

\seealso{
\code{\link{tvgarch}},
\code{\link[garchx]{garchx}},
\code{\link{tvgarchSim}}
}

\examples{
set.seed(12345)

## Simulate from a TV(1)-GARCH(1,1) model (default):
ySim <- tvgarchSim(n = 1500)

## Test of a TV(1)-GARCH(1,1) model:
yTest <- tvgarchTest(y = ySim)
orderG1 <- summary(yTest)

## Estimate a TV(1)-GARCH(1,1) model:
yEst <- tvgarch(y = ySim, order.g = orderG1)
}

\keyword{Time Varying Parameter Models}
\keyword{Nonlinear Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
