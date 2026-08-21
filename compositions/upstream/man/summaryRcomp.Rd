\name{summary.rcomp}
\alias{summary.rcomp}
\title{Summary of compositions in real geometry}
\description{
  Compute a summary of a composition based on real geometry.
}
\usage{
          \method{summary}{rcomp}( object, \dots ,robust=NULL)
          }
\arguments{
  \item{object}{an \code{\link{rcomp}} dataset of compositions }
  \item{\dots}{further arguments to \code{summary}}
  \item{robust}{A robustness description. See
    \link{robustnessInCompositions} for details. The option is currently
  not supported. If support is added the default will change to
  getOption(robust). }
}
\value{
  A matrix containing summary statistics.
  The value is the same as for the classical summary
  \code{\link[base]{summary}} applied to  a closed dataset.
}
\details{
  The data is applied a \code{\link{clo}} operation before the computation.
  Note that the statistics obtained will not keep any consistency 
  if computed with all the parts available or only with a subcomposition.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{rcomp}}, \code{\link{summary.aplus}}, \code{\link{summary.acomp}}
}
\examples{
data(SimulatedAmounts)
summary(rcomp(sa.lognormals))
summary(rcomp(sa.tnormals))

}
\keyword{multivariate}
