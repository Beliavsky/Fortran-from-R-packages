\name{ilrBase}
\alias{ilrBase}
%\alias{gsi.ilrBase}
\alias{ilrBaseList}
\title{The canonical basis in the clr plane used for ilr and ipt transforms.}
\description{
  Compute the basis of a clr-plane, to use with isometric log-ratio or planar transform of a (dataset of)
  compositions.
}
\usage{
          ilrBase( x=NULL , z=NULL , D = NULL, method = "basic" )
%         gsi.ilrBase(D)
          }
\arguments{
  \item{x}{ optional dataset or vector of compositions}
  \item{z}{ optional dataset or vector containing ilr or ipt coordinates}
  \item{D}{ number of parts of the simplex}
  \item{method}{ method to build the basis, one of "basic", "balanced", "optimal"
           "PBhclust", "PBmaxvar" or "PBangprox"}
}
\value{
  All methods give a matrix containing by columns the basis elements for the
  canonical basis of the clr-plane used for the ilr and ipt transform. Only one of the
  arguments \code{x}, \code{z} or \code{D} is needed
  to determine the dimension of the simplex.
  
  If you provide transformed data \code{z}, the function attempts to extract the basis 
  information from it with \code{\link{gsi.getV}}. Otherwise, the default compatible 
  ilr base matrix is created.
}
\details{
  Method "basic" computes a triangular Helmert matrix (corresponding to
  the original ilr transformation defined by Egozcue et al, 2003).
  In this case, \code{ilrBase} is a wrapper catching
  the answers of \code{gsi.ilrBase} and is to be
  used as the more convenient function.

  Method "balanced" returns an ilr matrix associated with a balanced partition,
  splitting the parts in groups as equal as possible. Transforms \code{ilr} and \code{ipt} computed
  with this basis are less affected by any component (as happens with "basic").

  The following methods are all data-driven and will fail if \code{x} is not given. 
  Some of these methods are extended to non-acomp datasets via the \code{cpt}
  general functionality. Use with care with non-acomp objects!

  Method "optimal" is a wrapper to \code{gsi.optimalilrBase}, providing the ilr basis
  with less influence of missing values. It is computed as a hierarchical
  cluster of variables, with parts previously transformed to
  1 (if the value is lost) or 0 (if it is recorded).
  
  Methods "PBhclust", "PBmaxvar" and "PBangprox" are principal balance methods (i.e.
  balances approximating principal components in different ways). These are all
  resolved by calls to \code{gsi.PrinBal}. Principal balances functionality should be
  considered beta!
}
\seealso{
  \code{\link{clr}},\code{\link{ilr}},\code{\link{ipt}}
}
\references{
  Egozcue J.J., V. Pawlowsky-Glahn, G. Mateu-Figueras and
  C. Barcel'o-Vidal (2003) Isometric logratio transformations for
  compositional data analysis. \emph{Mathematical Geology}, \bold{35}(3)
  279-300\cr

  \url{https://ima.udg.edu/Activitats/CoDaWork03/}
  }
\examples{
ilr(c(1,2,3))
ilrBase(D=2)
ilrBase(c(1,2,3))
ilrBase(z= ilr(c(1,2,3)) )
round(ilrBase(D=7),digits= 3)
ilrBase(D=7,method="basic")
ilrBase(D=7,method="balanced")
}
\keyword{multivariate}
