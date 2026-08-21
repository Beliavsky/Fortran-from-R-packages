\name{simplemissingplot}
\alias{simpleMissingSubplot}
\title{Ternary diagrams}
\description{
  Displaying compositions in ternary diagrams
}
\usage{
simpleMissingSubplot(loc, portions, labels=NULL,
          col=c("white","yellow","red","green","blue"), 
          ..., border="gray60", vertical=NULL, xpd=NA)
}
\arguments{
  \item{loc}{a vector of the form c(x1,x2,y1,y2) giving the drawing
    rectangle for the subplot in coordinates as in par("usr"). I.e. if
    the plot is logrithmic the base 10 logarithm is to be used:
    \describe{
      \item{x1}{left boundary of drawing rectangle}
      \item{y1}{lower boundary of drawing rectangle}
      \item{x2}{right boundary of drawing rectangle}
      \item{y2}{upper boundary of drawing rectangle}
    }
  }
  \item{portions}{The portions of different missing categories}
  \item{labels}{The labels for the categories. }
  \item{col}{The colors to plot the different categories.}
  \item{\dots}{further graphical parameters passed to
    \code{\link{text}}}
  \item{border}{The color to draw the borders of the rectangles.}
  \item{vertical}{Should a horizontal or a vertical plot be produced. If
  NULL the choice is done automatically according to the size of the
  recangle provided.}
  \item{xpd}{extended plot region. See \code{par("xpd")}.}
}
\details{
  This function is typically not called directly,however it could in
  principle be used to add to plots.
  The user will modify the function call only to modify the appearance
  of the missing plot.
  
  The labels are only
    plotted for nonzero portions. In this way it is always possible to
    realize the presence of a given missing type, even if it is a too small
  portion to be actually displayed. In case of overplotting of different
labels a further investigation using missingSummary should be used.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{
  \code{\link{plot.aplus}} }

\examples{
data(SimulatedAmounts)
plot(acomp(sa.missings))
plot(acomp(sa.missings),mp=~simpleMissingSubplot(c(0,0.1,0.2,1),
   missingInfo[c(1,3:5,2)], 
   c("Not Missing",paste("Missing Only:",cn),"Totally Missing"),
   col=c("gray","red","green","blue","darkgray"))
)
ms <- missingSummary(sa.missings)
for( i in 1:3 )
  simpleMissingSubplot(c(0.9+0.03*(i-1),0.9+0.03*i,0.2,1), ms[i,])
}
\keyword{hplot}
