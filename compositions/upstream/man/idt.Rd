\name{idt}
\alias{idt}
\alias{idt.default}
\alias{idt.acomp}
\alias{idt.ccomp}
\alias{idt.rcomp}
\alias{idt.aplus}
\alias{idt.rplus}
\alias{idt.rmult}
\alias{idt.factor}
\alias{idt.data.frame}
\alias{idtInv}
\alias{idtInv.default}
\alias{idtInv.acomp}
\alias{idtInv.ccomp}
\alias{idtInv.rcomp}
\alias{idtInv.aplus}
\alias{idtInv.rplus}
\alias{idtInv.rmult}
\alias{idtInv.factor}
\alias{idtInv.data.frame}
\title{Isometric default transform}
\description{
  Compute the isometric default transform of a vector (or dataset) of
  compositions or amounts in the selected class.
}
\usage{
          idt(x,...)
          \method{idt}{default}( x,... )
          \method{idt}{acomp}( x ,...)
          \method{idt}{rcomp}( x ,...)
          \method{idt}{aplus}( x ,...)
          \method{idt}{rplus}( x ,...)
          \method{idt}{rmult}( x ,...)
          \method{idt}{ccomp}( x ,...)
          \method{idt}{factor}( x ,...)
          \method{idt}{data.frame}( x ,...)
          idtInv(x,orig=gsi.orig(x),...)
          \method{idtInv}{default}( x ,orig=gsi.orig(x),...)
          \method{idtInv}{acomp}( x ,orig=gsi.orig(x), V=gsi.getV(x),...)
          \method{idtInv}{rcomp}( x ,orig=gsi.orig(x), V=gsi.getV(x),...)
          \method{idtInv}{aplus}( x ,orig=gsi.orig(x),...)
          \method{idtInv}{rplus}( x ,orig=gsi.orig(x),...)
          \method{idtInv}{ccomp}( x ,orig=gsi.orig(x),...)
          \method{idtInv}{rmult}( x ,orig=gsi.orig(x),...)
          \method{idtInv}{factor}( x ,orig=gsi.orig(x), V=gsi.getV(x),...)
          \method{idtInv}{data.frame}( x , orig=gsi.orig(x), ...)

          }
\arguments{
  \item{x}{a classed amount or composition, to be transformed with its
    isometric default transform, or its inverse; in case of the method for \code{\link{data.frame}}
    objects, the function attempts to track information about a previous class (in an attribute
    \code{origClass}, and if found, a transformation is tried with it; for factors, idt expands
    the factor according to the contrasts represented by \code{V}, or vice-versa.)}
  \item{\dots}{generic arguments past to underlying functions}
  \item{orig}{a compositional object which should be mimicked
             by the inverse transformation. It is the generic
	     argument. Typically the \code{orig} argument is stored as an attribute 
	     in \code{x} and will be extracted automatically by this function; if this
	     fails, \code{orig} can be set equal to the dataset that
	     was transformed in the first place.}
  \item{V}{matrix of (\emph{transposed, inverted}) logcontrasts; 
  together with \code{orig}, it defines the back-transformation.
   Typically the \code{V} argument is stored as an attribute 
	     in \code{x} and will be extracted automatically by this function; if this
	     fails, \code{V} must be manually set to the matrix V used in the idt/ilr/ipt
	     calculations. Argument not used in amounts or counts geometries.}
}
\value{
  A corresponding matrix of row-vectors containing the transforms. (Exception: idt.data.frame can return a data.frame if the input has no "origClass"-attribute)
}
\details{
  The general idea of this package is to analyse the same data with
  different geometric concepts, in a fashion as similar as possible. For each of the
  four concepts there exists an isometric transform expressing the geometry
  in a full-rank euclidean vector space. Such a transformation is computed 
  by \code{idt}. For \code{\link{acomp}} the transform is \code{\link{ilr}}, for
  \code{\link{rcomp}} it is \code{\link{ipt}}, for 
  \code{\link{aplus}} it is \code{\link{ilt}}, and for
  \code{\link{rplus}} it is \code{\link{iit}}. Keep in mind that the
  transform does not keep the variable names, since there is no guaranteed 
  one-to-one relation between the original parts and each transformed
  variable.
  \cr
  The inverse \code{idtInv} is intended to allow for an "easy" and automatic 
  back-transformation, without intervention of the user. The argument \code{orig} 
  (the one determining the behaviour of \code{idtInv} as a generic function) 
  tells the function which back-transformation should be applied, and
  gives the column names of \code{orig} to the back-transformed
  values of \code{x}. Therefore, it is very conventient to give the original classed
  data set used in the analysis as \code{orig}.
}
\author{R. Tolosana-Delgado, K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{backtransform}}, \code{\link{cdt}}, \code{\link{ilr}}, \code{\link{ipt}},
  \code{\link{ilt}}, \code{\link{cdtInv}}, \code{\link{ilrInv}}, \code{\link{iptInv}}, 
            \code{\link{iltInv}}, \code{\link{iitInv}}
}
\references{
    van den Boogaart, K.G. and R. Tolosana-Delgado (2008) "compositions": a unified 
  R package to analyze Compositional Data, \emph{Computers &
    Geosciences}, 34 (4), pages 320-338, \doi{10.1016/j.cageo.2006.11.017}.
}
\examples{
\dontrun{
# the idt is defined by
idt         <- function(x) UseMethod("idt",x)
idt.default <- function(x) x
idt.acomp   <- function(x) ilr(x) 
idt.rcomp   <- function(x) ipt(x) 
idt.aplus   <- ilt 
idt.rplus   <- iit 
}
idt(acomp(1:5))
idt(rcomp(1:5))
  data(Hydrochem)
  x = Hydrochem[,c("Na","K","Mg","Ca")]
  y = acomp(x)
  z = idt(y)
  y2 = idtInv(z,y)
  par(mfrow=c(2,2))
  for(i in 1:4){plot(y[,i],y2[,i])}
}
\keyword{multivariate}
