\name{variation}
\alias{variation}
\alias{variation.default}
\alias{variation.acomp}
\alias{variation.rcomp}
\alias{variation.aplus}
\alias{variation.rplus}
\alias{variation.rmult}
\alias{is.variation}
\title{Variation matrices of amounts and compositions}
\description{
  Compute the variation  matrix in the various approaches of compositional
 and amount  data analysis. Pay attention that this is not computing the variance or
 covariance matrix!
}
\usage{
    variation(x,\dots)
          \method{variation}{acomp}(x, \dots,robust=getOption("robust"))
          \method{variation}{rcomp}(x, \dots,robust=getOption("robust"))
          \method{variation}{aplus}(x, \dots,robust=getOption("robust"))
          \method{variation}{rplus}(x, \dots,robust=getOption("robust"))
          \method{variation}{rmult}(x, \dots,robust=getOption("robust"))
          is.variation(M, tol=1e-10)
          }
\arguments{
  \item{x}{a dataset, eventually of amounts or compositions}
  \item{\dots}{currently unused}
  \item{robust}{
    A description of a robust estimator. FALSE for the classical
    estimators. See \link{robustnessInCompositions} for
    further details.
  }
  \item{M}{a matrix, to check if it is a valid variation}
  \item{tol}{tolerance for the check}
}
\value{
  The variation matrix of x.
  
  For \code{is.variation}, a boolean saying if the matrix satisfies the conditions to be a variation matrix.
}
\details{
   The variation matrix was defined in the \code{\link{acomp}} context of
   analysis of compositions as the matrix of variances of all 
   possible log-ratios among components (Aitchison, 1986). The 
   generalization to rcomp objects is simply to reproduce the 
   variance of all possible differences between components. The 
   amount (\code{\link{aplus}}, \code{\link{rplus}}) and rmult objects 
   should not be treated with variation 
   matrices, because this was intended to skip the existence of a closure
   (which does not exist in the case of amounts).
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{cdt}}, \code{\link{clrvar2ilr}}, \code{\link{clo}},
\code{\link{mean.acomp}}, \code{\link{acomp}}, \code{\link{rcomp}},
\code{\link{aplus}}, \code{\link{rplus}} }

\examples{
data(SimulatedAmounts)
meanCol(sa.lognormals)
variation(acomp(sa.lognormals))
variation(rcomp(sa.lognormals))
variation(aplus(sa.lognormals))
variation(rplus(sa.lognormals))
variation(rmult(sa.lognormals))

}
\keyword{multivariate}
