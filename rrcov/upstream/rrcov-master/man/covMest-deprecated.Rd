\name{covMest}
\alias{covMest}
%
\title{Constrained M-Estimates of Location and Scatter}
\concept{High breakdown point}
\description{
    Computes constrained M-Estimates of multivariate location and scatter
    based on the translated biweight function (\sQuote{t-biweight}) using
    a High breakdown point initial estimate. The default initial estimate
    is the Minimum Volume Ellipsoid computed
    with \code{\link{CovMve}}. The raw (not reweighted) estimates are taken
    and the covariance matrix is standardized to determinant 1.
}
\usage{
covMest(x, cor=FALSE, r = 0.45, arp = 0.05, eps=1e-3,
    maxiter=120, control, t0, S0)
}
\arguments{
  \item{x}{a matrix or data frame. }
  \item{cor}{should the returned result include a correlation matrix?
    Default is \code{cor = FALSE}}.
  \item{r}{ required breakdown point.  Allowed values are between
    \code{(n - p)/(2 * n)} and 1 and the default is 0.45}
  \item{arp}{ asympthotic rejection point, i.e. the fraction of points
    receiving zero weight (see Rocke (1996)).  Default is \code{0.05}.}
  \item{eps}{ a numeric value specifying the relative precision of the solution of
    the M-estimate. Defaults to \code{1e-3}}
  \item{maxiter}{ maximum number of iterations allowed in the computation of
    the M-estimate. Defaults to 120}
  \item{control}{ a list with estimation options - same as these provided in the
    fucntion specification. If the control object is supplied, the parameters from it
    will be used. If parameters are passed also in the invocation statement, they will
    override the corresponding elements of the control object.}
  \item{t0}{ optional initial high breakdown point estimates of the location.
    If not supplied MVE will be used. }
  \item{S0}{ optional initial high breakdown point estimates of the
    scatter. If not supplied MVE will be used. }
}
\details{
  Rocke (1996) has shown that the S-estimates of multivariate location and scatter
  in high dimensions can be sensitive to outliers even if the breakdown point
  is set to be near 0.5. To mitigate this problem he proposed to utilize
  the translated  biweight (or t-biweight) method with a
  standardization step consisting of equating the median of \code{rho(d)}
  with the median under normality. This is then not an S-estimate, but is
  instead a constrained M-estimate. In order to make the smooth estimators
  to work, a reasonable starting point is necessary, which will lead reliably to a
  good solution of the estimator. In \code{covMest} the MVE computed by
  \code{\link{CovMve}} is used, but the user has the possibility to give her own
  initial estimates.

}
\note{
    The psi, rho and weight functions for the M estimation are encapsulated in a
    virtual S4 class \code{PsiFun} from which a \code{PsiBwt} class, implementing
    the translated biweight (t-biweight), is dervied. The base class  \code{PsiFun}
    contains also the M-iteration itself. Although not documented and not
    accessibale directly by the user these classes will form the bases for adding
    other functions (biweight, LWS, etc.) as well as S-estimates.
}
\author{Valentin Todorov \email{valentin.todorov@chello.at},

    (some code from C. Becker -
    http://www.sfb475.uni-dortmund.de/dienst/de/content/struk-d/bereicha-d/tpa1softw-d.html)
}
\value{
  An object of class \code{"mest"} which is basically a \code{\link{list}} with
  the following components. This class is "derived" from \code{"mcd"} so that
  the same generic functions -  \code{print},  \code{plot},  \code{summary} - can
  be used.
  NOTE: this is going to change - in one of the next revisions \code{covMest}
  will return an S4 class \code{"mest"} which is derived (i.e. \code{contains})
  form class \code{"cov"}.

  \item{center}{the final estimate of location.}
  \item{cov}{the final estimate of scatter.}
  \item{cor}{the estimate of the correlation matrix (only if
    \code{cor = TRUE}).}
  \item{mah}{mahalanobis distances of the observations using the
    M-estimate of the location and scatter.}
  \item{X}{the input data as a matrix.}
  \item{n.obs}{total number of observations.}
  \item{method}{character string naming the method (M-Estimates).}
  \item{call}{the call used (see \code{\link{match.call}}).}
}
\references{
    D.L.Woodruff and D.M.Rocke (1994) Computable robust estimation of multivariate location
    and shape on high dimension using compound estimators, \emph{Journal of the American
    Statistical Association}, \bold{89}, 888--896.

    D.M.Rocke (1996) Robustness properties of S-estimates of multivariate location and
    shape in high dimension, \emph{Annals of Statistics}, \bold{24}, 1327-1345.

    D.M.Rocke and D.L.Woodruff (1996) Identification of outliers in multivariate data \emph{Journal
    of the American Statistical Association}, \bold{91}, 1047--1061.

  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\seealso{
  \code{\link[robustbase]{covMcd}}
}
\keyword{robust}
\keyword{multivariate}
