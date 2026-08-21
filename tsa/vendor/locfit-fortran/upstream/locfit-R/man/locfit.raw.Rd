\name{locfit.raw}
\alias{locfit.raw}
\title{
Local Regression, Likelihood and Density Estimation.
}
\usage{
locfit.raw(x, y, weights=1, cens=0, base=0,
  scale=FALSE, alpha=0.7, deg=2, kern="tricube", kt="sph",
    acri="none", basis=list(NULL),
  deriv=numeric(0), dc=FALSE,
  family, link="default",
  xlim, renorm=FALSE,
  ev=rbox(),
  maxk=100, itype="default", mint=20, maxit=20, debug=0,
  geth=FALSE, sty="none")
}
\description{
  \code{locfit.raw} is an interface to Locfit using numeric vectors
  (for a model-formula based interface, use \code{\link{locfit}}).
  Although this function has a large number of arguments, most users
  are likely to need only a small subset.

  The first set of arguments (\code{x}, \code{y}, \code{weights},
  \code{cens}, and \code{base}) specify the regression
  variables and associated quantities.

  Another set (\code{scale}, \code{alpha}, \code{deg}, \code{kern},
    \code{kt}, \code{acri} and \code{basis}) control the amount of smoothing:
    bandwidth, smoothing weights and the local model. Most of these arguments
    are deprecated - they'll currently still work, but should be provided through
    the \code{lp()} model term instead.

  \code{deriv} and \code{dc} relate to derivative (or local slope)
    estimation.

  \code{family} and \code{link} specify the likelihood family.

  \code{xlim} and \code{renorm} may be used in density estimation.

  \code{ev} specifies the evaluation structure or set of evaluation points.

  \code{maxk},  \code{itype}, \code{mint}, \code{maxit} and \code{debug}
    control the Locfit algorithms, and will be rarely used.

  \code{geth} and \code{sty} are used by other functions calling
    \code{locfit.raw}, and should not be used directly.
}
\arguments{
\item{x}{
Vector (or matrix) of the independent variable(s). Can be constructed using the
\code{lp()} function.
}
\item{y}{
Response variable for regression models. For density families,
\code{y} can be omitted.
}
\item{weights}{
Prior weights for observations (reciprocal of variance, or sample size).  
}
\item{cens}{
Censoring indicators for hazard rate or censored regression. The coding
is \code{1} (or \code{TRUE}) for a censored observation, and
\code{0} (or \code{FALSE}) for uncensored observations.
}
\item{base}{
Baseline parameter estimate. If provided, the local regression model is
fitted as \eqn{Y_i = b_i + m(x_i) + \epsilon_i}, with Locfit estimating
the \eqn{m(x)} term. For regression models, this effectively subtracts
\eqn{b_i} from \eqn{Y_i}. The advantage of the \code{base} formulation
is that it extends to likelihood regression models.
}

\item{scale}{ Deprecated - see \code{\link{lp}()}. }
\item{alpha}{Deprecated - see \code{\link{lp}()}.
A single number (e.g. \code{alpha=0.7})
is interpreted as a nearest neighbor fraction. With two
componentes (e.g. \code{alpha=c(0.7,1.2)}), the first component
is a nearest neighbor fraction, and the second component is
a fixed component. A third component is the penalty term in locally
adaptive smoothing.
}
\item{deg}{
  Degree of local polynomial. Deprecated - see \code{\link{lp}()}.
}
\item{kern}{
  Weight function, default = \code{"tcub"}.
  Other choices are \code{"rect"}, \code{"trwt"}, \code{"tria"},
  \code{"epan"}, \code{"bisq"} and \code{"gauss"}. Choices may be restricted
  when derivatives are required; e.g. for confidence bands and some
  bandwidth selectors.
}
\item{kt}{
  Kernel type, \code{"sph"} (default); \code{"prod"}.
  In multivariate problems, \code{"prod"} uses a 
  simplified product model which speeds up computations.
}
\item{acri}{Deprecated - see \code{\link{lp}().}}
\item{basis}{User-specified basis functions.} 
%See \code{\link{lfbas}} for more details on this argument.}

\item{deriv}{
  Derivative estimation. If \code{deriv=1}, the returned fit will be
  estimating the derivative (or more correctly, an estimate of the
  local slope). If \code{deriv=c(1,1)} the second order derivative
  is estimated. \code{deriv=2} is for the partial derivative, with
  respect to the second variable, in multivariate settings.  }
\item{dc}{ Derivative adjustment.  }

\item{family}{
  Local likelihood family; \code{"gaussian"}; 
  \code{"binomial"}; \code{"poisson"}; \code{"gamma"} and \code{"geom"}.
  Density and rate estimation families are \code{"dens"}, \code{"rate"} and 
  \code{"hazard"} (hazard rate). If the family is preceded by a \code{'q'}
  (for example, \code{family="qbinomial"}), quasi-likelihood variance
  estimates are used. Otherwise, the residual variance (\code{\link{rv}})
  is fixed at 1. The default family is \code{"qgauss"} if a response
  \code{y} is provided; \code{"density"} if no response is provided.
}
\item{link}{
  Link function for local likelihood fitting. Depending on the family, 
  choices may be \code{"ident"}, \code{"log"}, \code{"logit"},
  \code{"inverse"}, \code{"sqrt"} and \code{"arcsin"}.
}

\item{xlim}{
  For density estimation, Locfit allows the density to be supported on
  a bounded interval (or rectangle, in more than one dimension).
  The format should be \code{c(ll,ul)} where \code{ll} is a vector of
  the lower bounds and \code{ur} the upper bounds. Bounds such as
  \eqn{[0,\infty)} are not supported, but can be effectively
  implemented by specifying a very large upper bound.
}
\item{renorm}{Local likelihood density estimates may not integrate
  exactly to 1. If \code{renorm=T}, the integral will be estimated
  numerically and the estimate rescaled. Presently this is implemented
  only in one dimension.
}

\item{ev}{
  The evaluation structure,
  \code{\link{rbox}()} for tree structures;
  \code{\link{lfgrid}()} for grids;
  \code{\link{dat}()} for data points;
  \code{\link{none}()} for none.
  A vector or matrix of evaluation points can also be provided,
  although in this case you may prefer to use the
  \code{\link{smooth.lf}()} interface to Locfit.
  Note that arguments \code{flim}, \code{mg} and \code{cut} are now
  given as arguments to the evaluation structure function, rather
  than to \code{locfit.raw()} directly (change effective 12/2001).
}

\item{maxk}{
  Controls space assignment for evaluation structures.
  For the adaptive evaluation structures, it is impossible to be sure
  in advance how many vertices will be generated. If you get
  warnings about `Insufficient vertex space', Locfit's default assigment
  can be increased by increasing \code{maxk}. The default is \code{maxk=100}.
}
\item{itype}{
  Integration type for density estimation. Available methods include
  \code{"prod"}, \code{"mult"} and \code{"mlin"}; and \code{"haz"} for
  hazard rate estimation problems. The available integration methods
  depend on model specification (e.g. dimension, degree of fit). By
  default, the best available method is used.
}
\item{mint}{
  Points for numerical integration rules. Default 20.
}
\item{maxit}{
  Maximum iterations for local likelihood estimation. Default 20.
}
\item{debug}{If > 0; prints out some debugging information.}

\item{geth}{Don't use!}
\item{sty}{ Deprecated - see \code{\link{lp}()}. }
}
\value{
An object with class "locfit". A standard set of methods for printing, ploting, etc. these objects is provided.
}
\references{
Loader, C., (1999) Local Regression and Likelihood.
}
\keyword{smooth}
