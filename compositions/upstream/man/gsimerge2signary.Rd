\name{gsi.merge2signary}
\alias{gsi.merge2signary}
\alias{gsi.ilrBase2signary}
\alias{gsi.optimalilrBase}
\alias{gsi.buildilrBase}
\alias{gsi.signary2ilrBase}
\alias{gsi.OrderIlr}
\title{Auxiliary functions to compute user-defined ilr and ipt transforms.}
\description{
  Compute the basis of a clr-plane, to use with isometric log-ratio or planar transform of a (dataset of)
  compositions.
}
\usage{
          gsi.merge2signary( M )
          gsi.ilrBase2signary( V )
          gsi.optimalilrBase( x )
          gsi.buildilrBase( W=c(1,-1) )
          gsi.signary2ilrBase( W=c(1,-1) )
          gsi.OrderIlr( V )
          }
\arguments{
  \item{M}{ a merge structure (as explained in \code{hclust}) }
  \item{x}{ a compositional data set }
  \item{W}{ a signary matrix (as explained below) defining a partition }
  \item{V}{ a matrix of change of basis from clr/cpt to ilr/ipt }
}
\value{
  These functions will not be usually called on themselves, but 
  through their wrappers, mainly \code{ilrBase}. Functions
  \code{gsi.merge2signary} and \code{gsi.ilrBase2signary} return 
   a signary matrix (as explained in "details"), \code{gsi.optimalilrBase} 
  returns a merge structure (as epxlained in \code{hclust}), and
  \code{gsi.buildilrBase} (and its alias \code{gsi.signary2ilrBase}) 
  returns an ilr matrix. These functions are thought to be called 
  sequentially. \cr
  Apart, \code{gsi.OrderIlr} reorders both parts and coordinates to 
  ease dendrogram-like representations (see \code{\link{CoDaDendrogram}}).
}
\details{
  A signary matrix is a matrix with the same shape as an ilr matrix, 
  but containing only +1, 0 or -1 values (thus, it is a kind of "extended 
  binary"). If the value W[i,j]= +1, then part "i" is involved in 
  coordinate "j" in the numerator; if W[i,j]=-1, it is involed
  in the denominator, and if W[i,j]=0 then part "i" does not take 
  part in coordinate "j".

  Functions \code{gsi.merge2signary} and \code{gsi.buildilrBase} are 
  intended to compute \code{ilrBase} matrices associated to user-defined 
  partition structures. Function \code{gsi.ilrBase2signary} offers the 
  inverse functionality.

  Function \code{gsi.OrderIlr} returns a list with two elements: 
  "ilrBase" and "order". The first one contains the ilr basis with 
  coordinates reordered in decreasing number of involved parts (so, 
  all parts are involved in the first coordinate, and only two in the 
  last). The second one gives a permutation of the parts so that 
  involved parts in each coordinate are always together. 
  Note that ilrBase does not have its parts permuted!
}
\note{It is better not to use gsi.* functions directly since they are internal
  functions of the package. Use their wrappers.}
\author{Raimon Tolosana-Delgado, K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{ilrBase}},\code{\link{ipt}},\code{\link{ilr}},
  \url{https://ima.udg.edu/Activitats/CoDaWork03/}
}
\references{
  Egozcue J.J., V. Pawlowsky-Glahn, G. Mateu-Figueras and
  C. Barcel'o-Vidal (2003) Isometric logratio transformations for
  compositional data analysis. \emph{Mathematical Geology}, \bold{35}(3)
  279-300\cr
  }
\examples{
m <- matrix(data=c(-1,-2,
                  -3,-4,
                   1, 2),ncol=2,nrow=3,byrow=TRUE)
w <- gsi.merge2signary(m)
w
V=gsi.buildilrBase(w)
V
gsi.ilrBase2signary(V)
gsi.OrderIlr(V)
}
\keyword{internal}
