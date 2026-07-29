\name{VaR_ES_bounds_analytical}
\alias{crude_VaR_bounds}
\alias{VaR_bounds_hom}
\alias{dual_bound}
\title{``Analytical'' Best and Worst Value-at-Risk for Given Marginals}
\description{
  Compute the best and worst Value-at-Risk (VaR) for given marginal
  distributions with an ``analytical'' method.
}
\usage{
## ``Analytical'' methods
crude_VaR_bounds(level, qF, d = NULL, ...)
VaR_bounds_hom(level, d, method = c("Wang", "Wang.Par", "dual"),
               interval = NULL, tol = NULL, ...)
dual_bound(s, d, pF, tol = .Machine$double.eps^0.25, ...)
}
\arguments{
  \item{level}{confidence level \eqn{\alpha}{alpha} for VaR and ES
    (e.g., 0.99).}
  \item{qF}{\code{d}-list containing the marginal quantile functions.
    In the homogeneous case, \code{qF} can also be a single function.}
  \item{d}{dimension (number of risk factors; \eqn{\ge 2}{ >= 2}). For
    \code{crude_VaR_bounds()}, \code{d} only needs to be given in
    the homogeneous case in which \code{qF} is a \code{\link{function}}.}
  \item{method}{\code{\link{character}} string. \code{method = "Wang"} and
    \code{method = "Wang.Par"}
    apply the approach of McNeil et al. (2015,
    Proposition 8.32) for computing best (i.e., smallest) and
    worst (i.e., largest) VaR. The
    latter method assumes Pareto margins and thus does
    not require numerical integration.
    \code{method = "dual"} applies the dual bound approach as
    in Embrechts et al. (2013, Proposition 4)
    for computing worst VaR (no value for the best
    VaR can be obtained with this approach and thus
    \code{\link{NA}} is returned for the best VaR).}
  \item{interval}{initial interval (a \code{\link{numeric}(2)}) for
    computing worst VaR. If not provided, these are the defaults chosen:
    \describe{
      \item{\code{method = "Wang"}:}{initial
	interval is \eqn{[0,(1-\alpha)/d]}{[0,(1-alpha)/d]}.}
      \item{\code{method = "Wang.Par"}:}{initial
	interval is \eqn{[c_l,c_u]}, where \eqn{c_l} and \eqn{c_u}
	are chosen as in Hofert et al. (2015).}
      \item{\code{method = "dual"}:}{in this case, no good defaults are known.
	Note that the lower endpoint of the initial interval has to be
	sufficiently large in order for the the inner root-finding algorithm
	to find a root; see Details.}
  }}
  \item{tol}{tolerance for
    \code{\link{uniroot}()}
    for computing worst VaR. This defaults (for \code{tol = NULL}) to
    \eqn{2.2204*10^{-16}}{2.2204*10^{-16}}
    for \code{method = "Wang"} or \code{method = "Wang.Par"} (where a
    smaller tolerance is crucial) and to \code{\link{uniroot}()}'s
    default \code{.Machine$double.eps^0.25} otherwise. Note that for
    \code{method = "dual"}, \code{tol} is used for both the outer
    and the inner root-finding procedure.}
  \item{s}{dual bound evaluation point.}
  \item{pF}{marginal loss distribution function (homogeneous case only).}
  \item{\dots}{
    \describe{
      \item{\code{crude_VaR_bounds()}:}{ellipsis argument passed
        to (all provided) quantile functions.}
      \item{\code{VaR_bounds_hom()}:}{case \code{method = "Wang"}
        requires the quantile function \code{qF()} to be provided
        and additional arguments passed via the ellipsis argument are passed
	on to
	% the underlying objective function Wang_h() from which they
	% are passed on to
	the underlying \code{\link{integrate}()}. For
        \code{method = "Wang.Par"}
	the ellipsis argument must contain the parameter \code{shape}
	(the shape parameter \eqn{\theta>0}{theta>0} of the Pareto distribution).
	For \code{method = "dual"}, the ellipsis argument must contain the
	distribution function \code{pF()} and the initial interval
	\code{interval} for the outer root finding procedure (not for
	\code{d = 2}); additional arguments are passed on to the underlying
	\code{\link{integrate}()} for computing the dual bound
	\eqn{D(s)}.}
      \item{\code{dual_bound()}:}{ellipsis argument is passed to the
	underlying \code{\link{integrate}()}.}
    }}
}
\value{
  \code{crude_VaR_bounds()} returns crude lower and upper bounds for
  VaR at confidence level \eqn{\alpha}{alpha} for any
  \eqn{d}-dimensional model with marginal quantile functions
  specified by \code{qF}.

  \code{VaR_bounds_hom()} returns the best and worst VaR at
  confidence level \eqn{\alpha}{alpha} for \eqn{d} risks with equal
  distribution function specified by the ellipsis \code{...}.

  \code{dual_bound()} returns the value of the dual bound \eqn{D(s)} as
  given in Embrechts, Puccetti, \enc{Rüschendorf}{Rueschendorf}
  (2013, Eq. (12)).
}
\details{
  For \code{d = 2}, \code{VaR_bounds_hom()} uses the method of
  Embrechts et al. (2013,
  Proposition 2). For \code{method = "Wang"} and \code{method = "Wang.Par"}
  the method presented in McNeil et al. (2015, Prop. 8.32) is
  implemented; this goes back to Embrechts et al. (2014, Prop. 3.1; note that
  the published version of this paper contains typos for both bounds).
  This requires
  one \code{\link{uniroot}()} and, for the generic \code{method = "Wang"},
  one \code{\link{integrate}()}. The critical part for the
  generic \code{method = "Wang"} is the lower endpoint of the initial
  interval for \code{\link{uniroot}()}. If the (marginal)
  distribution function has finite first moment, this can be taken as
  0. However, if it has infinite first moment, the lower endpoint has to
  be positive (but must lie below the unknown root). Note that the upper
  endpoint \eqn{(1-\alpha)/d}{(1-alpha)/d} also happens to be a
  root and thus one needs a proper initional interval containing the
  root and being stricticly contained in
  \eqn{(0,(1-\alpha)/d}{(1-alpha)/d)}.
  In the case of Pareto margins, Hofert et al. (2015) have
  derived such an initial (which is used by
  \code{method = "Wang.Par"}).
  Also note that the chosen smaller default tolerances for
  \code{\link{uniroot}()} in case of \code{method = "Wang"} and
  \code{method = "Wang.Par"} are crucial for obtaining reliable
  VaR values; see Hofert et al. (2015).

  For \code{method = "dual"} for computing worst VaR, the method
  presented of Embrechts et al. (2013, Proposition 4) is implemented.
  This requires two (nested) \code{\link{uniroot}()}, and an
  \code{\link{integrate}()}. For the inner root-finding procedure to
  find a root, the lower endpoint of the provided initial
  \code{interval} has to be \dQuote{sufficiently large}.

  Note that these approaches for computing the
  VaR bounds in the homogeneous case are numerically non-trivial;
  see the source code and \code{vignette("VaR_bounds",
    package = "qrmtools")}
  for more details. As a
  rule of thumb, use \code{method = "Wang"} if you have to (i.e., if the
  margins are not Pareto) and \code{method = "Wang.Par"} if you can (i.e.,
  if the margins are Pareto). It is not recommended to use
  (the numerically even more challenging) \code{method = "dual"}.
}
\author{Marius Hofert}
\references{
  Embrechts, P., Puccetti, G., \enc{Rüschendorf}{Rueschendorf}, L.,
  Wang, R. and Beleraj, A. (2014).  An Academic Response to Basel
  3.5. \emph{Risks} \bold{2}(1), 25--48.

  Embrechts, P., Puccetti, G. and \enc{Rüschendorf}{Rueschendorf}, L. (2013).
  Model uncertainty and VaR aggregation. \emph{Journal of Banking &
    Finance} \bold{37}, 2750--2764.

  McNeil, A. J., Frey, R. and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.

  Hofert, M., Memartoluie, A., Saunders, D. and Wirjanto, T. (2017).
  Improved Algorithms for Computing Worst
  Value-at-Risk. \emph{Statistics & Risk Modeling}
  or, for an earlier version, \url{https://arxiv.org/abs/1505.02281}.
}
\seealso{
  \code{\link{RA}()}, \code{\link{ARA}()}, \code{\link{ABRA}()}
  for empirical solutions in the inhomogeneous case.

  \code{vignette("VaR_bounds", package = "qrmtools")}
  for more example calls, numerical challenges
  encoutered and a comparison of the different methods for computing
  the worst (i.e., largest) Value-at-Risk.
}
\examples{
## See ?rearrange
}
\keyword{programming}