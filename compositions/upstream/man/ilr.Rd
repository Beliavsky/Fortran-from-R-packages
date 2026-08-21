\name{ilr}
\alias{ilr}
\alias{ilrInv}
\title{Isometric log ratio transform}
\description{
  Compute the isometric log ratio transform of a (dataset of)
  composition(s), and its inverse.
}
\usage{
          ilr( x , V = ilrBase(x) ,...)
          ilrInv( z , V = ilrBase(z=z),..., orig=gsi.orig(z))
          }
\arguments{
  \item{x}{a composition, not necessarily closed}
  \item{z}{the ilr-transform of a composition}
  \item{V}{a matrix, with columns giving the chosen basis of the
    clr-plane}
  \item{\dots}{generic arguments. not used.}
  \item{orig}{a compositional object which should be mimicked 
             by the inverse transformation. It is especially used to
	     reconstruct the names of the parts.}
}
\value{
  \code{ilr} gives the isometric log ratio transform,
  \code{ilrInv} gives closed compositions with the given ilr-transforms
}
\details{
  The ilr-transform maps a composition in the D-part Aitchison-simplex
  isometrically to a D-1 dimensonal euclidian vector. The data can then
  be analysed in this transformation by all classical multivariate
  analysis tools. However the interpretation of the results may be
  difficult, since there is no one-to-one relation between the original parts 
  and the transformed variables.\cr
  
  The isometric logratio transform is given by
  \deqn{ ilr(x) := V^t clr(x)  }
  with \code{\link{clr}}(x) the centred log ratio transform and
  \eqn{V\in R^{d \times (d-1)}}{V} a matrix which columns form an orthonormal 
  basis of the clr-plane. A default matrix \eqn{V} is given by
  \code{ilrBase(\var{D})}.
  
   
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}
\seealso{
  \code{\link{clr}},\code{\link{alr}},\code{\link{apt}}, \code{\link{ilrBase}}
}
\references{
  Egozcue J.J., V. Pawlowsky-Glahn, G. Mateu-Figueras and
  C. Barcel'o-Vidal (2003) Isometric logratio transformations for
  compositional data analysis. \emph{Mathematical Geology}, \bold{35}(3)
  279-300\cr
  Aitchison, J, C. Barcel'o-Vidal, J.J. Egozcue, V. Pawlowsky-Glahn
  (2002) A consise guide to the algebraic geometric structure of the
  simplex, the sample space for compositional data analysis, \emph{Terra
    Nostra}, Schriften der Alfred Wegener-Stiftung, 03/2003
  \cr
  \url{https://ima.udg.edu/Activitats/CoDaWork03/}
  }
\examples{
(tmp <- ilr(c(1,2,3)))
ilrInv(tmp)
ilrInv(tmp) - clo(c(1,2,3)) # 0
data(Hydrochem)
cdata <- Hydrochem[,6:19]
pairs(ilr(cdata))
ilrBase(D=3) 
}
\keyword{multivariate}
