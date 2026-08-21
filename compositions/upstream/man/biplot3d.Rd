\name{biplot3D}
\alias{biplot3D}
\alias{biplot3D.default}
\alias{biplot3D.princomp}
\title{Three-dimensional biplots, based on package rgl}
\description{
  Plots variables and cases in the same plot, based on a principal
  component analysis. 
}
\usage{
biplot3D(x,...)
\method{biplot3D}{default}(x,y,var.axes=TRUE,col=c("green","red"),cex=c(2,2),
            xlabs = NULL, ylabs = NULL, expand = 1,arrow.len = 0.1,
            ...,add=FALSE)
 \method{biplot3D}{princomp}(x,choices=1:3,scale=1,...,
            comp.col=1,comp.labs=paste("Comp.",1:3),
            scale.scores=lambda[choices]^(1-scale),
            scale.var=scale.comp, scale.comp=sqrt(lambda[choices]), 
            scale.disp=1/scale.comp)
}
\arguments{
  \item{x}{princomp object or matrix of point locations to be drawn
    (typically, cases)}
  \item{choices}{Which principal components should be used?}
  \item{scale}{a scaling parameter like in \code{\link{biplot}}}
  \item{scale.scores}{a vector giving the scaling applied to the scores}
  \item{scale.var}{a vector giving the scaling applied to the variables}
  \item{scale.comp}{a vector giving the scaling applied to the unit length
    of each component}
  \item{scale.disp}{a vector giving the scaling of the display in the
    directions of the components}
  \item{comp.col}{color to draw the axes of the components, defaults to black}
  \item{comp.labs}{labels for the components}
  \item{\dots}{further plotting parameters as defined in %\code{\link[rgl]{rgl.material}}
  \code{rgl::material3d}}
  \item{y}{matrix of  second point/arrow-head locations (typically, variables)}
  \item{var.axes}{logical, TRUE draws arrows and FALSE points for y}
  \item{col}{vector/list of two elements the first giving the
    color/colors for the first data set and the second giving
    color/colors for the second data set.}
  \item{cex}{vector/list of two elements the first giving the
    size for the first data set and the second giving
    size for the second data set.}
  \item{xlabs}{labels to be plotted at x-locations}
  \item{ylabs}{labels to be plotted at y-locations}
  \item{expand}{the relative expansion of the y data set with respect to
    x}
  \item{arrow.len}{The length of the arrows as defined in \code{\link{arrows3D}}}
  \item{add}{logical, adding to existing plot or making a new one?}
}
\value{
  the 3D plotting coordinates of the tips of the arrows of the variables 
  displayed, returned invisibly
}
\details{
  This "biplot" is a triplot, relating data, variables and principal
  components. The relative scaling of the components is still
  experimental, meant to mimic the behavior of classical
  biplots.
}

\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{gsi}}
}
\examples{
data(SimulatedAmounts)
pc <- princomp(acomp(sa.lognormals5))
pc
summary(pc)
plot(pc)      #plot(pc,type="screeplot")
biplot(pc)
if(requireNamespace("rgl", quietly = TRUE)) {
  biplot3D(pc)
} ## this function requires package 'rgl'

}
\keyword{multivariate}
