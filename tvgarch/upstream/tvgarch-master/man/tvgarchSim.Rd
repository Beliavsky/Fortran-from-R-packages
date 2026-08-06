\name{tvgarchSim}
\alias{tvgarchSim}

\title{Simulate from a univariate TV-GARCH-X model}

\description{
Simulate from a univariate multiplicative TV(s)-GARCH(p,q,r)-X model.
}

\usage{
tvgarchSim(n, order.g = 1, order.h = c(1,1,0),
           intercept.g = 1.2, size = 5, speed = 25, location = 0.5, xtv = NULL, 
           intercept.h = 0.2, arch = 0.1, garch = 0.8, asym = NULL, xreg = NULL,
           opt = 0, as.zoo = TRUE, verbose = FALSE, innovations = NULL)
}

\arguments{
  \item{n}{\code{integer}.}
  \item{order.g}{\code{integer} \code{vector} of length s indicating the number of locations in each transition function of the TV component.}
  \item{order.h}{\code{integer} \code{vector} of the form c(p,q,r). The first entry controls the GARCH order, the second the ARCH order and the third the asymmetry order of the GARCH-X component.}
  \item{intercept.g}{\code{NULL} or numeric with the value of the intercept in the TV component.}
  \item{size}{\code{NULL} or numeric \code{vector} with the values of the size coefficients.}
  \item{speed}{\code{NULL} or numeric \code{vector} with the values of the \code{speed} coefficients.}
  \item{location}{\code{NULL} or numeric \code{vector} with the values of the location coefficients.}
  \item{xtv}{\code{NULL} or numeric \code{vector}, time series or zoo object to include as the transition variable in the TV component. If \code{NULL}, calendar time, scaled between 0 and 1, is used as the transition variable.}
  \item{opt}{\code{integer} indicating whether the \code{speed} parameter in the TV component should be scaled. If 0, no scaling; if 1, \code{speed}/sd(\code{xtv}); if 2, exp(\code{speed}).}
  \item{intercept.h}{numeric with the value of the intercept in the GARCH-X component.}
  \item{arch}{\code{NULL} or numeric \code{vector} with the values of the ARCH-coefficients.}
  \item{garch}{\code{NULL} or numeric \code{vector} with the values of the GARCH-coefficients.}
  \item{asym}{\code{NULL} or numeric \code{vector} with the values of the asymmetry-coefficients.}
  \item{xreg}{\code{NULL} or numeric \code{vector} with the values of the X-term.}
  \item{as.zoo}{\code{logical}. If \code{TRUE}, then the returned result is of class \code{\link[zoo]{zoo}}.}
  \item{verbose}{\code{logical}, if \code{TRUE}, the conditional variance and innovations are also returned.}
  \item{innovations}{\code{NULL} or numeric \code{vector} with the innovations. If \code{NULL}, then standard normal innovations are generated with \code{\link{rnorm}}.}
}

\value{An object of class 'zoo' (if \code{as.zoo = TRUE}), otherwise a vector or a matrix 
(if \code{verbose = TRUE}), with the simulated values.}

\author{Susana Campos-Martins}

\seealso{
\code{\link{tvgarch}},
\code{\link[garchx]{garchx}},
\code{\link[zoo]{zoo}}
}

\examples{
set.seed(123)

## Simulate from a TV(1)-GARCH(1,1) model (default):
ySim1 <- tvgarchSim(n = 1500)

## Simulate from a TV(2)-GARCH(1,1) model:
ySim2 <- tvgarchSim(n = 1500, order.g = c(1,2), size = c(0.5,-0.4), 
                    speed = c(1.5,2), location = c(0.2, 0.5,0.8))

## Simulate from a GARCH(1,1) model:
ySim3 <- tvgarchSim(n = 1500, order.g = NULL)

## Simulate from a TV(1)-GARCH(1,1,1)-X model:
ySim4 <- tvgarchSim(n = 1500, order.h = c(1,1,1), asym = 0.025, xreg = ySim3^2)
}

\keyword{Time Varying Parameter Models}
\keyword{Nonlinear Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
