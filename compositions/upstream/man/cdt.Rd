\name{cdt}
\alias{cdt}
\alias{cdt.default}
\alias{cdt.acomp}
\alias{cdt.rcomp}
\alias{cdt.aplus}
\alias{cdt.rplus}
\alias{cdt.rmult}
\alias{cdt.ccomp}
\alias{cdt.factor}
\alias{cdt.data.frame}
\alias{cdtInv}
\alias{cdtInv.default}
\alias{cdtInv.acomp}
\alias{cdtInv.rcomp}
\alias{cdtInv.aplus}
\alias{cdtInv.rplus}
\alias{cdtInv.rmult}
\alias{cdtInv.ccomp}
\alias{cdtInv.factor}
\alias{cdtInv.data.frame}
\title{Centered default transform}
\description{
  Compute the centered default transform of a (data set of)
  compositions or amounts (or its inverse).
}
\usage{
          cdt(x,...)
          \method{cdt}{default}( x ,...)
          \method{cdt}{acomp}( x ,...)
          \method{cdt}{rcomp}( x ,...)
          \method{cdt}{aplus}( x ,...)
          \method{cdt}{rplus}( x ,...)
          \method{cdt}{rmult}( x ,...)
          \method{cdt}{ccomp}( x ,...)
          \method{cdt}{factor}( x ,...)
          \method{cdt}{data.frame}( x ,...)
          cdtInv(x,orig=gsi.orig(x),...)
          \method{cdtInv}{default}( x ,orig=gsi.orig(x),...)
          \method{cdtInv}{acomp}( x ,orig=gsi.orig(x),...)
          \method{cdtInv}{rcomp}( x ,orig=gsi.orig(x),...)
          \method{cdtInv}{aplus}( x ,orig=gsi.orig(x),...)
          \method{cdtInv}{rplus}( x ,orig=gsi.orig(x),...)
          \method{cdtInv}{rmult}( x ,orig=gsi.orig(x),...)
          \method{cdtInv}{ccomp}( x ,orig=gsi.orig(x),...)
          \method{cdtInv}{factor}( x ,orig=gsi.orig(x),...)
          \method{cdtInv}{data.frame}( x ,orig=gsi.orig(x),...)
          }
\arguments{
  \item{x}{a classed (matrix of) amount or composition, to be transformed with its
    centered default transform, or its inverse; in case of the method for \code{\link{data.frame}}
    objects, the function attempts to track information about a previous class (in an attribute
    \code{origClass}, and if found, a cdt/Inv transformation is tried with it; for factors, 
    cdt expands the factor in indicator values for each category, or vice-versa.)}
  \item{\dots}{generic arguments past to underlying functions.}
  \item{orig}{a compositional object which should be mimicked
             by the inverse transformation. It is used to determine the
	     backtransform to be used and eventually to
	     reconstruct the names of the parts. It is the generic
	     argument. Typically this argument is extracted from \code{x},
	     but if this fails you can give the data set that
	     has be transformed in the first place.}
}
\value{
  A corresponding matrix or vector containing the transforms. (Exception: cdt.data.frame can return
  a data.frame if the input has no "origClass"-attribute)
}
\details{
  The general idea of this package is to analyse the same data with
  different geometric concepts, in a fashion as similar as possible. For each of the
  four concepts there exists a unique transform expressing the geometry
  in a linear subspace, keeping the relation to the variables. This
  unique transformation is computed by \code{cdt}. For
  \code{\link{acomp}} the transform is \code{\link{clr}}, for
\code{\link{rcomp}} it is \code{\link{cpt}}, for 
  \code{\link{aplus}} it is \code{\link{ilt}}, and for
  \code{\link{rplus}} it is \code{\link{iit}}. Each component of the result
  is identified with a unit vector in the direction of the corresponding
  component of the original composition or amount. Keep in mind that the
  transform is not necessarily surjective and thus variances in the
  image space might be singular.

  
}
\author{R. Tolosana-Delgado, K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{backtransform}}, \code{\link{idt}}, \code{\link{clr}}, \code{\link{cpt}}, \code{\link{ilt}}, \code{\link{iit}}
}
\references{
    van den Boogaart, K.G. and R. Tolosana-Delgado (2008) "compositions": a unified 
  R package to analyze Compositional Data, \emph{Computers &
    Geosciences}, 34 (4), pages 320-338, \doi{10.1016/j.cageo.2006.11.017}.
}
\examples{
\dontrun{
# the cdt is defined by
cdt         <- function(x) UseMethod("cdt",x)
cdt.default <- function(x) x
cdt.acomp   <- clr 
cdt.rcomp   <- cpt 
cdt.aplus   <- ilt 
cdt.rplus   <- iit 
}
x <- acomp(1:5)
(ds <- cdt(x))
cdtInv(ds,x)
(ds <- cdt(rcomp(1:5)))
cdtInv(ds,rcomp(x))
  data(Hydrochem)
  x = Hydrochem[,c("Na","K","Mg","Ca")]
  y = acomp(x)
  z = cdt(y)
  y2 = cdtInv(z,y)
  par(mfrow=c(2,2))
  for(i in 1:4){plot(y[,i],y2[,i])}

}
\keyword{multivariate}
