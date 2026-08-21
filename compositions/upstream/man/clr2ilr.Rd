\name{clr2ilr}
\alias{clr2ilr}
\alias{ilr2clr}
\alias{clrvar2ilr}
\alias{ilrvar2clr}
\alias{clrvar2variation}
\alias{variation2clrvar}
\alias{is.clrvar}
\alias{is.ilrvar}
\title{Convert between clr and ilr, and between cpt and ipt. }
\description{
  Compute the centered log ratio transform of a (dataset of) from
  isometric log-ratio transform(s) and its inverse. Equivalently, 
  compute centered and isometric planar transforms from each other.
  Acts in vectors and in bilinear forms. For bilinear forms, 
  transform between variation-form from clr-form.
}
\usage{
clr2ilr( x , V=ilrBase(x=x) )
ilr2clr( z , V=ilrBase(z=z), x=gsi.orig(z) )
clrvar2ilr( varx , V=ilrBase(D=ncol(varx)) )
ilrvar2clr( varz , V=ilrBase(D=ncol(varz)+1) ,x=NULL)
clrvar2variation(Sigma)
variation2clrvar(TT)
is.clrvar(M, tol=1e-10)
is.ilrvar(M, tol=1e-10)
          }
\arguments{
  \item{x}{the clr/cpt-transform of composition(s) (in the ilr2-routines
    provided only to give column names.)}
  \item{z}{the ilr/ipt-transform of composition(s)}
  \item{varx, Sigma}{variance or covariance matrix of clr/cpt-transformed
    compositions}
  \item{varz}{variance or covariance matrix of ilr/ipt-transformed
    compositions}
  \item{V}{a matrix with columns giving the chosen  basis of the clr-plane}
  \item{TT}{variation matrix}
  \item{M}{a matrix, to check if it is a valid variance}
  \item{tol}{tolerance for the check}

}
\value{
  \code{clr2ilr} gives the ilr/ipt transform of the same composition(s),\cr
  \code{ilr2clr} gives the clr/cpt transform of the same
  composition(s),\cr
  \code{clrvar2ilr} gives the variance-/covariance-matrix of the ilr/ipt transform of the same compositional data set,\cr
  \code{ilrvar2clr} and \code{clrvar2variation} give the variance-/covariance-matrix of the clr/cpt
  transform of the same compositional data set.\cr
  \code{variation2clrvar} gives the variation matrix from the clr-covariance matrix\cr
  \code{is.*var} check if the given matrix satisfies the conditions to be an ilr-variance 
  resp. a clr-variance
}
\details{
  These functions perform a matrix multiplication with \code{V} in an
  appropriate way. 
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{
  \code{\link{variation}}, \code{\link{ilr}}, \code{\link{ipt}},
  \code{\link{clr}}, \code{\link{cpt}}
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
  

}
\examples{
data(SimulatedAmounts)
ilrInv(clr2ilr(clr(sa.lognormals)))-clo(sa.lognormals)
clrInv(ilr2clr(ilr(sa.lognormals)))-clo(sa.lognormals)
ilrvar2clr(var(ilr(sa.lognormals)))-var(clr(sa.lognormals))
clrvar2ilr(var(cpt(sa.lognormals)))-var(ipt(sa.lognormals))
variation(acomp(sa.lognormals))
clrvar2variation(var(acomp(sa.lognormals)))
}
\keyword{multivariate}
