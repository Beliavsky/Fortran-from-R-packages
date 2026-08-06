\name{CovMve}
\alias{CovMve}
%
\title{ Robust Location and Scatter Estimation via MVE }
\description{
     Computes a robust multivariate location and scatter estimate with a high
     breakdown point, using the \sQuote{MVE} (Minimum Volume Ellipsoid)
     estimator.
}
\usage{
CovMve(x, alpha = 1/2, nsamp = 500, seed = NULL, trace = FALSE, control)
}
\arguments{
  \item{x}{a matrix or data frame. }
  \item{alpha}{numeric parameter controlling the size of the subsets
    over which the determinant is minimized, i.e., \code{alpha*n}
    observations are used for computing the determinant.  Allowed values
    are between 0.5 and 1 and the default is 0.5.}
  \item{nsamp}{ number of subsets used for initial estimates or \code{"best"}
    or \code{"exact"}.  Default is \code{nsamp = 500}.  For
    \code{nsamp="best"} exhaustive enumeration is done, as long as the
    number of trials does not exceed 5000. For \code{"exact"},
    exhaustive enumeration will be attempted however many samples are
    needed.  In this case a warning message will be displayed saying
    that the computation can take a very long time.}
  \item{seed}{starting value for random generator. Default is \code{seed = NULL}}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
  \item{control}{ a control object (S4) of class \code{\link{CovControlMve-class}}
    containing estimation options - same as these provided in the fucntion
    specification. If the control object is supplied, the parameters from it
    will be used. If parameters are passed also in the invocation statement, they will
    override the corresponding elements of the control object.}
}
\details{
  This function computes the minimum volume ellipsoid estimator
  of location and scatter and returns an S4 object of class
  \code{\link{CovMve-class}} containing the estimates.

  The approximate estimate is
  based on a subset of size \code{alpha*n} with an enclosing ellipsoid of smallest volume.
  The mean of the best found subset provides the raw estimate of the location,
  and the rescaled covariance matrix is the raw estimate of scatter. The rescaling of
  the raw covariance matrix is by \code{median(dist)/qchisq(0.5, p)} and this scale factor
  is returned in the slot \code{raw.cnp2}. Currently no finite sample corrction factor is applied.
  The Mahalanobis distances of all observations from the location estimate for the
  raw covariance matrix are calculated, and those points within the 97.5% point
  under Gaussian assumptions are declared to be good. The final (reweightd) estimates are the
  mean and rescaled covariance of the good points. The reweighted covariance matrix is
  rescaled by \code{1/pgamma(qchisq(alpha, p)/2, p/2 + 1)/alpha}
  (see Croux and Haesbroeck, 1999) and this scale factor is returned
  in the slot \code{cnp2}.

  The search for the approximate solution is made over ellipsoids determined by the
  covariance matrix of \code{p+1} of the data points and applying
  a simple but effective improvement of the subsampling procedure
  as described in Maronna et al. (2006), p. 198.
  Although there exists no formal proof of this improvement (as for MCD and LTS),
  simulations show that it can be recommended as an approximation of the MVE.
}
\note{
    Main reason for implementing the MVE estimate was that it is the recommended
    initial estimate for S estimation (see Maronna et al. (2006), p. 199) and will
    be used by default in \code{\link{CovMest}} (after removing the correction
    factors from the covariance matrix and rescaling to determinant 1).
}
\value{
  An S4 object of class \code{\link{CovMve-class}} which is a subclass of the
  virtual class \code{\link{CovRobust-class}}.
}
\references{
  P. J. Rousseeuw and A. M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection.} Wiley.

  C. Croux and G. Haesbroeck (1999).
    Influence function and efficiency of the minimum covariance determinant scatter matrix estimator.
    \emph{Journal of Multivariate Analysis},
    \bold{71}, 161--190.

    R. A. Maronna, D. Martin and V. Yohai (2006). \emph{Robust Statistics: Theory and Methods}.
    Wiley, New York.

  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} and
Matias Salibian-Barrera \email{matias@stat.ubc.ca}
}
\seealso{
  \code{\link[MASS]{cov.rob}} from package \pkg{MASS}
}
\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
CovMve(hbk.x)

## the following three statements are equivalent
c1 <- CovMve(hbk.x, alpha = 0.75)
c2 <- CovMve(hbk.x, control = CovControlMve(alpha = 0.75))
## direct specification overrides control one:
c3 <- CovMve(hbk.x, alpha = 0.75,
             control = CovControlMve(alpha=0.95))
c1
}
\keyword{robust}
\keyword{multivariate}
