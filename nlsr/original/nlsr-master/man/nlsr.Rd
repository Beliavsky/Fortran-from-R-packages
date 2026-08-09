\name{nlsr}
\alias{nlsr}
\title{nlsr function}
\usage{
  nlsr(formula = NULL,  data = NULL,  start = NULL,  control = NULL,
    trace = FALSE, subset = NULL, lower = -Inf, upper = Inf,  weights = NULL,
    ...)
}
\arguments{
\item{formula}{The modeling formula. Looks like 'y~b1/(1+b2*exp(-b3*T))'}

\item{data}{a data frame containing data for variables
used in the formula that are NOT the parameters. This
data may also be defined in the parent frame i.e.,
'global' to this function}

\item{start}{MUST be a named vector with all elements present
e.g., start=c(b1=200, b2=50, b3=0.3)}

\item{control}{a list of control parameters. See nlsr.control().}

\item{trace}{TRUE for console output during execution (default FALSE)}

\item{subset}{an optional vector specifying a subset of observations
to be used in the fitting process. NOT used currently by nlxb()
or nlfb() and will throw an error if present and not NULL.}

\item{lower}{a vector of lower bounds on the parameters.
If a single number, this will be applied to all parameters
Default \code{-Inf}.}

\item{upper}{a vector of upper bounds on the parameters. If a single number,
this will be applied to all parameters. Default \code{Inf}.}

\item{weights}{A vector of fixed weights. The objective function that will be
minimized is the sum of squares where each residual is multiplied by the
square root of the corresponding weight. Default \code{NULL} implies
unit weights.}

\item{...}{additional data needed to evaluate the modeling functions}
}
\value{
A solution object of type \code{nls}
}
\description{
Provides class nls solution to a nonlinear least squares solution
using the Nash Marquardt tools.
}
