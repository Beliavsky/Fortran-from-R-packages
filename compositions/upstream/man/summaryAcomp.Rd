\name{summary.acomp}
\alias{summary.acomp}
\title{Summarizing a compositional dataset in terms of ratios}
\description{
  Summaries in terms of compositions are quite different from classical
  ones. Instead of analysing each variable individually, we must
  analyse each pair-wise ratio in a log geometry.
}
\usage{
          \method{summary}{acomp}( object, \dots ,robust=getOption("robust"))
          }
\arguments{
  \item{object}{a data matrix of compositions, not necessarily closed}
  \item{\dots}{not used, only here for generics}
  \item{robust}{A robustness description. See
    \link{robustnessInCompositions} for details. The parameter can be
    null for avoiding any estimation.}
}
\value{
  The result is an object of type \code{"summary.acomp"}
  \item{mean}{the \code{\link{mean.acomp}} composition}
  \item{mean.ratio}{a matrix containing the geometric mean of the
    pairwise ratios}
  \item{variation}{the variation matrix of the dataset (\code{{\link{variation.acomp}}})}
  \item{expsd}{a matrix containing the  one-sigma factor for
    each ratio, computed as \code{exp(sqrt(variation.acomp(W)))}. To
    obtain a two-sigma-factor, one has to take its squared value (power 1.96, actually).}
  \item{invexpsd}{the inverse of the preceding one, giving the reverse bound.
    Additionally, it can be "almost" intepreted as a correlation coefficient, 
    with values near one indicating high proportionality between the components.}
  \item{min}{a matrix containing the minimum of each of the pairwise ratios}
  \item{q1}{a matrix containing the 1-Quartile of each of the pairwise ratios}
  \item{median}{a matrix containing the median of each of the pairwise ratios}
  \item{q1}{a matrix containing the 3-Quartile of each of the pairwise ratios}
  \item{max}{a matrix containing the maximum of each of the pairwise ratios}

}
\details{
  It is quite difficult to summarize a composition in a consistent and
  interpretable way. We tried to provide such a summary here, based on the
  idea of the variation matrix.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, R. Tolosana-Delgado}
\seealso{
  \code{\link{acomp}}
}
\references{
  Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.
}
\examples{
data(SimulatedAmounts)
summary(acomp(sa.lognormals))

}
\keyword{multivariate}
