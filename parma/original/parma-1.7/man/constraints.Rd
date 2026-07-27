\name{constraints}
\alias{ineqfun.turnover.min} 
\alias{ineqjac.turnover.min}
\alias{ineqfun.bsturnover.min}
\alias{ineqjac.bsturnover.min}
\alias{ineqfun.turnover.opt} 
\alias{ineqjac.turnover.opt}
\alias{ineqfun.bsturnover.opt}
\alias{ineqjac.bsturnover.opt}
\alias{ineqfun.variance.opt}
\alias{ineqjac.variance.opt}
\alias{ineqfun.variance.min}
\alias{ineqjac.variance.min}
\title{
NLP custom constraint functions
}
\description{
Provides a number of custom constraints and their jacobians for use with the NLP 
representation (both minimum risk and the fractional problem).
}
\usage{
ineqfun.turnover.min(w, optvars, uservars)
ineqjac.turnover.min(w, optvars, uservars)
ineqfun.bsturnover.min(w, optvars, uservars)
ineqjac.bsturnover.min(w, optvars, uservars)
ineqfun.turnover.opt(w, optvars, uservars)
ineqjac.turnover.opt(w, optvars, uservars)
ineqfun.bsturnover.opt(w, optvars, uservars)
ineqjac.bsturnover.opt(w, optvars, uservars)
ineqfun.variance.opt(w, optvars, uservars)
ineqjac.variance.opt(w, optvars, uservars)
ineqfun.variance.min(w, optvars, uservars)
ineqjac.variance.min(w, optvars, uservars)
}
\arguments{
  \item{w}{ The decision weight vector.}
  \item{optvars}{ Problem specific list (not for use by user).}
  \item{uservars}{ User specific list with some required parameters to be set 
  (see details).}
}
\details{
Functions preceded be \dQuote{ineqfun} denote the functions which may be
passed to the ineqfun option in the \code{parmaspec} function, while 
\dQuote{eqfun} to the eqfun option. Functions which include \dQuote{jac} in the
first part of the name denote the jacobians of the equivalent constraints and
should be passed to the ineqgrad or eqgrad options in \code{parmaspec}. Functions
ending with \dQuote{min} denote formulations for use with the minrisk type problems
while functions ending with \dQuote{opt} denote the fractional risk formulation.\cr
For the simple turnover constraint, the uservars list must contain an entry
called \sQuote{wold} denoting the previous vector of weights with which the
comparison will be made. Additionally, an entry called  \sQuote{turnover} 
is required which denotes the (positive) value for the maximum turnover.\cr
For the buy and sell turnover constraint (bsturnover), there should instead
be (positive) \sQuote{buyturnover} and \sQuote{sellturnover} entries in the
uservars list. Finally note that when using this type of constraints in a fractional 
programming setup, care should be taken that the combination of bounds, turnover 
limits and the forecast return vector do not result in a negative expected return 
in which case the problem is not solvable.\cr
The variance constraint allows the targeting of a maximum acceptable variance.
The extra arguments which must be passed to the uservars list are \sQuote{Cov} 
for the asset covariance matrix and \sQuote{varbound} representing the maximum 
acceptable upper variance.
}
\value{
Used internally by the NLP solver. The fun return a scalar, while the jac return
a matrix with n.cols equal to the length of the decision vector (which may be
greater than the length of the weights as in the fractional problem which 
contains the fractional multiplier as well as other formulations which have 
additional decision variables).
}
\author{
Alexios Galanos
}
\keyword{optimize}
