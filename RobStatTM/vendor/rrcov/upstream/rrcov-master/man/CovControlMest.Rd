\name{CovControlMest}
\alias{CovControlMest}
\title{ Constructor function for objects of class "CovControlMest" }
\concept{High breakdown point}
\description{
    This function will create a control object \code{CovControlMest}
    containing the control parameters for \code{CovMest}
}
\usage{
CovControlMest(r = 0.45, arp = 0.05, eps = 0.001, maxiter = 120)
}
\arguments{
    \item{r}{a numeric value specifying the required
    breakdown point.  Allowed values are between
    \code{(n - p)/(2 * n)} and 1 and the default is 0.45}
    \item{arp}{a numeric value specifying the asympthotic
    rejection point, i.e. the fraction of points receiving zero
    weight (see Rocke (1996)).  Default is \code{0.05} }
    \item{eps}{a numeric value specifying the
    relative precision of the solution of the M-estimate.
    Defaults to \code{1e-3}}
    \item{maxiter}{maximum number of iterations allowed
    in the computation of the M-estimate. Defaults to 120 }
}
%\details{}
\value{
A \code{CovControlMest} object
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
%\note{}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
%\seealso{}
\examples{
    ## the following two statements are equivalent
    ctrl1 <- new("CovControlMest", r=0.4)
    ctrl2 <- CovControlMest(r=0.4)

    data(hbk)
    CovMest(hbk, control=ctrl1)
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
