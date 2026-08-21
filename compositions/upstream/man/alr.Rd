\name{alr}
\alias{alr}
\alias{alrInv}
\title{Additive log ratio transform}
\description{
  Compute the additive log ratio transform of a (dataset of)
  composition(s), and its inverse.
}
\usage{
          alr( x ,ivar=ncol(x), ... )
          alrInv( z, ...,orig=gsi.orig(z))
          }
\arguments{
  \item{x}{a composition, not necessarily closed}
  \item{z}{the alr-transform of a composition, thus a (D-1)-dimensional
    real vector}
  \item{\dots}{generic arguments. not used.}
  \item{orig}{a compositional object which should be mimicked 
             by the inverse transformation. It is especially used to
	     reconstruct the names of the parts.}
  \item{ivar}{The column to be used as denominator variable. Unfortunately
    not yet supported in alrInv. The default works even if x is a vector.}
}
\value{
  \code{alr} gives the additive log ratio transform; accepts a compositional dataset
  \code{alrInv} gives a closed composition with the given alr-transform; accepts a dataset
}
\details{
  The alr-transform maps a composition in the D-part Aitchison-simplex
  non-isometrically to a D-1 dimensonal euclidian vector, treating the
  last part as common denominator of the others. The data can then
  be analysed in this transformation by all classical multivariate
  analysis tools not relying on a distance. The interpretation of
  the results is relatively simple, since the relation to the original D-1
  first parts is preserved. However distance is an extremely relevant
  concept in most types of analysis, where a \code{\link{clr}} or
  \code{\link{ilr}} transformation should be preferred.\cr
  
  The additive logratio transform is given by
  \deqn{ alr(x)_i := \ln\frac{x_i}{x_D}  }{alr(x)_i := \ln(x_i/x_D)}.
   
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{clr}},\code{\link{ilr}},\code{\link{apt}},
  \url{https://ima.udg.edu/Activitats/CoDaWork03/}
}
\references{
  Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.
}
\examples{
(tmp <- alr(c(1,2,3)))
alrInv(tmp)
unclass(alrInv(tmp)) - clo(c(1,2,3)) # 0
data(Hydrochem)
cdata <- Hydrochem[,6:19]
pairs(alr(cdata),pch=".")
}
\keyword{multivariate}
