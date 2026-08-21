\name{ipt}
\alias{ipt}
\alias{iptInv}
\alias{uciptInv}
\title{Isometric planar transform}
\description{
  Compute the isometric planar  transform of a (dataset of)
  composition(s) and its inverse.
}
\usage{
          ipt( x , V = ilrBase(x),... )
          iptInv( z , V = ilrBase(z=z),...,orig=gsi.orig(z))
          uciptInv( z , V = ilrBase(z=z),...,orig=gsi.orig(z) )
          }
\arguments{
  \item{x}{a composition or a data matrix of compositions, not necessarily closed}
  \item{z}{the ipt-transform of a composition or a data matrix of
    ipt-transforms of compositions}
  \item{V}{a matrix with columns giving the chosen basis of the clr-plane}
  \item{\dots}{generic arguments. not used.}
  \item{orig}{a compositional object which should be mimicked 
             by the inverse transformation. It is especially used to
	     reconstruct the names of the parts.}

}
\value{
  \code{ipt} gives the centered planar transform,
  \code{iptInv} gives closed compositions with with the given ipt-transforms,
  \code{uciptInv} unconstrained iptInv does the same as iptInv but
  sets illegal values to NA rather than giving an error. This is a
  workaround to allow procedures not honoring the constraints of the
  space.  
  }
\details{
  The ipt-transform maps a composition in the D-part real-simplex
  isometrically to a D-1 dimensonal euclidian vector. Although the
  transformation does not reach the whole \eqn{R^{D-1}}, resulting covariance
  matrices are typically of full rank.
  \cr 
  The data can then
  be analysed in this transformation by all classical multivariate
  analysis tools. However, interpretation of results may be
  difficult, since the
  transform does not keep the variable names, given that there is no
  one-to-one relation between the original parts and each transformed variables. See
  \code{\link{cpt}} and \code{\link{apt}} for alternatives. \cr
    
  The isometric planar transform is given by
    \deqn{ ipt(x) := V^t cpt(x)  }
  with \code{\link{cpt}}(x) the centred planar transform and
  \eqn{V\in R^{d \times (d-1)}}{V} a matrix which columns form an orthonormal 
  basis of the clr-plane. A default matrix \eqn{V} is given by
  \code{ilrBase(\var{D})}

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{ilr}},\code{\link{ilrBase}}, \code{\link{cpt}}
}
\references{
    van den Boogaart, K.G. and R. Tolosana-Delgado (2008) "compositions": a unified 
  R package to analyze Compositional Data, \emph{Computers &
    Geosciences}, 34 (4), pages 320-338, \doi{10.1016/j.cageo.2006.11.017}.

}
\examples{
(tmp <- ipt(c(1,2,3)))
iptInv(tmp)
iptInv(tmp) - clo(c(1,2,3)) # 0
data(Hydrochem)
cdata <- Hydrochem[,6:19]
pairs(ipt(cdata)) 
}
\keyword{multivariate}
