\name{tvgarchObj}

\alias{tv}
\alias{tvObj}
\alias{garchObj}
\alias{tvgarchObj}

\title{Auxiliary functions}

\description{
Auxiliary functions used in the estimation of the univariate and multivariate 
TV(s)-GARCH(p,q,r)-X model. Not intended for the average user.
}

\usage{
  tv(speed, location, xtv = NULL, n = NULL, opt = 0, 
  order.g = NULL, as.zoo = TRUE, verbose = FALSE)
  tvObj(par.g, fixed.par.g, xtv, opt, order.g, fixed.h, y, iter0, flag)
  garchObj(par.h, xreg, order.h, fixed.g, y, flag)
  tvgarchObj(par, fixed.par.g, y, order.g, xtv, opt, iter.fit.h, flag)
}

\arguments{
  \item{speed}{\code{NULL} or numeric \code{vector} with the values of the 
  speed coefficients.}
  \item{location}{\code{NULL} or numeric \code{vector} with the values of the 
  location coefficients.}
  \item{xtv}{\code{NULL} or numeric \code{vector}, time series or 
  \code{\link[zoo]{zoo}} 
  object to include as the transition variable in the TV component. If 
  \code{NULL}, a continuous variable bounded between 0 and 1 for \code{n} 
  observations is constructed and used as the transition variable.}
  \item{n}{\code{integer} indicating the number of observations of the 
  continuous transition variable bounded between 0 and 1.}
  \item{opt}{ \code{integer}, indicates whether the speed parameter in the TV 
  component should be scaled. If 0, no scaling; if 1, speed/sd(xtv); if 2, 
  exp(speed). For function \code{tv()}, the default is 0.}
  \item{order.g}{a scalar in \code{tv()} and an \code{integer} \code{vector} 
  of length s in \code{tvObj()} indicating the number of locations in each 
  transition function of the TV component. For function \code{tv()}, it defaults 
  to \code{NULL}.}
  \item{as.zoo}{\code{logical}. If \code{TRUE}, then the returned result is of 
  class \code{\link[zoo]{zoo}}.}
  \item{verbose}{\code{logical}, if \code{TRUE}, the values of not only the 
  logistic transition function but also the transition variable are returned.}
  \item{par.g}{numeric \code{vector} with the values of the parameters in the TV 
  component. If \code{iter0}=TRUE, \code{par.g} takes the form 
  c(\code{intercept.g}, \code{size}, \code{speed}, \code{location}); if 
  \code{iter0}=\code{FALSE}, then \code{par.g}=c(\code{speed},\code{size},
  \code{location}) and the values of the fixed parameters are provided using 
  \code{fixed.par.g}).
  }
  \item{fixed.par.g}{\code{NULL} or numeric \code{vector} with the values of the 
  parameters fixed in the TV component, i.e., intercept.g.}
  \item{par}{\code{NULL} or numeric \code{vector} with the values of the 
  parameters in the TV-GARCH-X model.}
  \item{fixed.h}{numeric \code{vector}, time series or \code{\link[zoo]{zoo}} 
  containing the values of GARCH-X component).}
  \item{y}{numeric \code{vector}, time series or \code{\link[zoo]{zoo}} object.}
  \item{iter0}{\code{logical}. If \code{FALSE}, some parameters in the TV 
  component are fixed during the iterative estimation.}
  \item{flag}{\code{integer}. If 0, returns a numeric \code{vector} with the 
  values of the objective function; if 1 returns the the value of the objective 
  function; if 2, returns the fitted variance components.}
  \item{par.h}{numeric \code{vector} with the values of the parameters in the 
  GARCH-X component.}
  \item{order.h}{\code{integer} \code{vector} of the form c(p,q,r). The first 
  entry controls the GARCH order, the second the ARCH order and the third the 
  asymmetry order of the GARCH-X component.}
  \item{xreg}{numeric \code{vector}, times eries or zoo object to include as 
  covariates in the GARCH-X component.}
  \item{fixed.g}{numeric \code{vector}, time series or \code{\link[zoo]{zoo}} 
  containing the values of TV component).}
  \item{iter.fit.h}{a list of class 'garchx'.}
}

\value{
The values of the objective function or fitted variance components.
}

\author{Susana Campos-Martins}

\seealso{
\code{\link{tvgarch}}, 
\code{\link{fitted.tvgarch}}, 
\code{\link{residuals.tvgarch}} 
}

\keyword{Time Varying Parameter Models}
\keyword{Nonlinear Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
