\name{barplot.acomp}
\alias{barplot.acomp}
\alias{barplot.rcomp}
\alias{barplot.aplus}
\alias{barplot.rplus}
\alias{barplot.ccomp}
\title{Bar charts of amounts}
\description{
  Compositions and amounts dispalyed as bar plots. 
}
\usage{
\method{barplot}{acomp}(height,...,legend.text=TRUE,beside=FALSE,total=1,
plotMissings=TRUE,missingColor="red",missingPortion=0.01)
\method{barplot}{rcomp}(height,...,legend.text=TRUE,beside=FALSE,total=1,
plotMissings=TRUE,missingColor="red",missingPortion=0.01)
\method{barplot}{aplus}(height,...,legend.text=TRUE,beside=TRUE,total=NULL,
plotMissings=TRUE,missingColor="red",missingPortion=0.01)
\method{barplot}{rplus}(height,...,legend.text=TRUE,beside=TRUE,total=NULL,
plotMissings=TRUE,missingColor="red",missingPortion=0.01)
\method{barplot}{ccomp}(height,...,legend.text=TRUE,beside=FALSE,total=1,
plotMissings=TRUE,missingColor="red",missingPortion=0.01)
          }
\arguments{
  \item{height}{an acomp, rcomp, aplus, or rplus object giving amounts
    to be displayed}
  \item{\dots}{further graphical parameters as in
    \code{\link[graphics]{barplot}}}
  \item{legend.text}{same as legend.text in \code{\link[graphics]{barplot}}}
  \item{beside}{same as beside in \code{\link[graphics]{barplot}}}
  \item{total}{The total to be used in displaying the composition,
    typically 1, 100 or the number of parts. If NULL no normalisation
    takes place.}
  \item{plotMissings}{logical: shall missings be annotate in the plot}
  \item{missingColor}{color to draw missings}
  \item{missingPortion}{The space portion to be reserved for missings}
}
\value{A numeric vector (or matrix, when \code{beside = TRUE}) giving
  the coordinates of all the bar midpoints drawn, as in
  \code{\link[graphics]{barplot}} }
\details{
  These functions are essentially light-weighted wrappers for
  \code{\link[graphics]{barplot}}, just adding an adequate default
  behavior for each of the scales. The missingplot functionality will
  work well with the default settings.

  If \code{plotMissings} is true, there will be an additional portion
  introduced, which is not counted in the total. This might make the
  plots looking less nice, however they make clear to the viewer that it
  is by no means clear how the rest of the plot should be interpreted
  and that the missing value really casts some unsureness on the rest of
  the data.  
  
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{acomp}}, \code{\link{rcomp}},   \code{\link{rplus}}
  \code{\link{aplus}}, \code{\link{plot.acomp}},
  \code{\link{boxplot.acomp}}


}
\examples{
data(SimulatedAmounts)
barplot(mean(acomp(sa.lognormals[1:10,])))
barplot(mean(rcomp(sa.lognormals[1:10,])))
barplot(mean(aplus(sa.lognormals[1:10,])))
barplot(mean(rplus(sa.lognormals[1:10,])))

barplot(acomp(sa.lognormals[1:10,]))
barplot(rcomp(sa.lognormals[1:10,]))
barplot(aplus(sa.lognormals[1:10,]))
barplot(rplus(sa.lognormals[1:10,]))

barplot(acomp(sa.tnormals))
}
\keyword{multivariate}
