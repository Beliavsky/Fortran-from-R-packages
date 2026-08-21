\name{perturbe}
\alias{perturbe}
\alias{+.acomp}
\alias{-.acomp}
\title{Perturbation of compositions}
\description{
  The perturbation is the addition operation in the Aitchison geometry
  of the simplex.
}
\usage{
perturbe(x,y)
## Methods for class "acomp"
## x + y
## x - y
##   - x
}
\arguments{
  \item{x}{compositions of class \code{\link{acomp}}}
  \item{y}{compositions of class \code{\link{acomp}}}
}
\value{
  An \code{\link{acomp}} vector or matrix.  }
\details{
  The perturbation is the basic addition operation of the Aitichson
  simplex as a vector space. It is defined by:
  \deqn{ (x+y)_i = clo( (x_i  y_i)_i )_i }
  \code{perturbe} and \code{+} compute this operation. The only
  difference is that \code{+} checks the class of its argument, while
  \code{perturbe} does not check the type of the arguments and can thus
  directly be applied to a composition in any form (unclassed, acomp,
  rcomp).\cr
  The \code{-} operation is the inverse of the addition in the usual way
  and defined by:
  \deqn{(x-y)_i:=clo( (x_i/y_i)_i )_i}
  and as unary operation respectively as: 
  \deqn{(-x)_i:=clo( (1/y_i)_i )_i}
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}
\seealso{
  \code{\link{acomp}}, \code{\link{*.aplus}}, \code{\link{+.rplus}}
}
\references{
  Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.
}
\examples{
tmp <- -acomp(1:3)
tmp + acomp(1:3)


}
\keyword{multivariate}
