\name{gpdbiv}
\alias{gpdbiv}
\title{Implements Bivariate POT Method}
\description{
  Returns an object of class \code{"gpdbiv"} representing the fit
  of a bivariate POT (peaks over thresholds) model for joint
  excesses over thresholds. 
}
\usage{
gpdbiv(data1 = NA, data2 = NA, u1 = NA, u2 = NA, ne1 = NA, ne2 = NA,
    global = FALSE, method = "BFGS", \dots)
}
\arguments{
  \item{data1}{first data vector}
  \item{data2}{second data vector}
  \item{u1}{threshold for \code{data1} (either this or \code{ne1}
    must be given but not both)}
  \item{u2}{threshold for \code{data2} (either this or \code{ne2}
    must be given but not both)}
  \item{ne1}{number of upper extremes to be used for \code{data1}
    (either this or \code{u1} must be given but not both)}
  \item{ne2}{number of upper extremes to be used for \code{data2}
    (either this or \code{u2} must be given but not both)}
  \item{global}{should a global maximisation of the likelihood with
    respect to marginal and dependence parameters be undertaken. The
    default alternative is a two-stage local fit where first the
    marginal parameters are estimated and then the dependence
    parameter. This is much faster than a global fit.}
  \item{method}{the optimization method (see \code{optim}).
    The argument has been created (as distinct from \dots) in order
    to make the \code{"BFGS"} method the default, as the default
    used by \code{optim} is not recommended for the one-dimensional
    optimizations that occur when \code{global = FALSE}.}
  \item{\dots}{other arguments passed to \code{optim}}
}
\value{
  An object of class \code{"gpdbiv"} representing the fit and including
  parameter estimates and standard errors. 
}
\details{
  This function implements a model suggested by Richard Smith (see
  references below). The marginal excess distributions are GPD
  distributions, as suggested by univariate EVT and implemented in gpd.  
  The dependence specification is known as the logistic or Gumbel
  dependence structure, but it would be easy to program alternatives. 
}
\references{
  Multivariate Threshold Methods, Richard L. Smith,
  in \emph{Extreme Value Theory and Applications}, ed. J. Galambos,
  published by Kluwer, pages 225-248, 1994.
  
  Markov Chain Models for Threshold Exceedances,
  R.L. Smith, J.A. Tawn, S.G. Coles,
  \emph{Biometrika} \bold{84}, 249-268, 1997. 
}
\seealso{\code{\link{gpd}}, \code{\link{plot.gpdbiv}}, 
  \code{\link{interpret.gpdbiv}}}

\examples{
data(bmw) ; data(siemens)
out <- gpdbiv(-bmw, -siemens, ne1 = 100, ne2 = 100)
interpret.gpdbiv(out, 0.05, 0.05)
\dontrun{plot(out)} 
}
\keyword{models}
