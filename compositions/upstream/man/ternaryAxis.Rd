\name{ternaryAxis}
\alias{ternaryAxis}
\title{Axis for ternary diagrams}
\description{
  Displaying compositions in ternary diagrams
}
\usage{ternaryAxis(side=1:3,at=seq(0.2,0.8,by=0.2),
                        labels=if(is.list(at)) lapply(at,format) else format(at),
                        ...,
                        tick=TRUE,pos=0,
                        font.axis=par("font.axis"),
                        font.lab=par("font.lab"),
                        lty="solid",lwd=1,
                        len.tck=0.025,dist.lab=0.03,
                        dist.axis=0.03,
                        lty.tck="solid",
                        col.axis=par("col.axis"),
                        col.lab=par("col.lab"),
                        cex.axis=par("cex.axis"),
                        cex.lab=par("cex.lab"),
                        Xlab=NULL,Ylab=NULL,Zlab=NULL,small=TRUE,
                        xpd=NA,aspanel=FALSE)}
\arguments{
  \item{side}{a vector giving the sides to draw the axis on. 1=under the
    plot, 2=the upper right axis, 3=the upper left axis. -1 is the
    portion axis of the first component, -2 is the portion axis of the
    second component, -3 is the portion axis of the third component.
    An empty vector
    or 0 suppresses axis plotting, but still plots the Xlab, Ylab and Zlab
    parameters.}
  \item{at}{ a vector or a list of vectors giving the positions of the
    tickmarks. }
  \item{labels}{a vector giving the labels or a list of things that
    can serve as graphics annotations. Each element of the list is than
    sean as the labels for one of axes. IMPORTANT: if plotting formulae
    enclose the list of labels into a list.}
  \item{tick}{a logical whether to draw the tickmark lines}
  \item{pos}{the portion of the opposite component to draw the axis on.
    Proportion axss shrinks, when pos>0 !
  }
  \item{font.axis}{the font for the axis annotations}
  \item{font.lab}{the font for the variable labels}
  \item{lty}{the line type of the axis line. (see \code{\link{par}}). NA
    supresses plotting.}
  \item{lty.tck}{the line type of the tickmarks. NA suppresses plotting.}
  \item{len.tck}{the line length of the tickmarks.}
  \item{dist.axis}{the distance of the variable labels from the
    axes. Positve values point outward from the plot.}
  \item{dist.lab}{the distance of the axes labels from the
    axes. Positve values point outward from the plot.}
  \item{lwd}{the line widths of axis line and tickmarks.
    (see \code{\link{par}})}
  \item{col.axis}{the color to plot the axis line, the tickmarks and the
    axes labels.}
  \item{col.lab}{the color to plot the variable labels.}
  \item{cex.axis}{The character size to plot the axes labels. (see
\code{\link{par}})}
  \item{cex.lab}{The character size for the variable labels}
  \item{Xlab}{the label for the lower left component.}
  \item{Ylab}{the label for the lower right component.}
  \item{Zlab}{the label for the upper component.}
  \item{small}{wether to plot the lower labels under the corners}
  \item{xpd}{Extended plotting region. See (see \code{\link{par}}).}
  \item{aspanel}{Is this called as a slave to acutally plot the axis
    (TRUE), or as a user level function to instatiate the axis (FALSE). }
  \item{\dots}{further graphical that might be of use for other
    functions, but are silently ignored here}
}
\details{
  This function has two uses. If called with \code{aspanel=TRUE} it
  acutally draws the axes to a panel. In other cases it tries to modify
  the axes argument of the current plot to add the axis. I.e. it will
  force a replotting of the plot with the new axes settings. Thus an old
  axes is removed.

  To ensure that various axes can be drawn with various parameters most
  of the arguments can take a vector or list of the same length as
  \code{side}
  providing the different parameters for each of the axes to be drawn.

  There are two types of axes: Proportion axes (1:3) and portions axes
  (-1:-3). The best place to draw a Proportion axes is pos=0, which is
  the standard for axis in ternary diagrams. Portion axes are best drawn
  at \code{pos=0.5} in the middle of the plot. 
}
\author{K.Gerald v.d. Boyogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}
\seealso{
  \code{\link{plot.aplus}}, \code{\link{plot3D}} (for 3D plot),
   \code{\link{kingTetrahedron}} (for 3D-plot model export),
  \code{\link{qqnorm.acomp}},\code{\link{boxplot.acomp}} }

\examples{
data(SimulatedAmounts)
plot(acomp(sa.lognormals),axes=TRUE)
ternaryAxis(side=1:3,pos=0,col.axis="red",col.lab="green")
ternaryAxis(side=1:3,at=1:9/10,
            labels=expression(9:1,4:1,7:3,3:2,1:1,2:3,3:7,1:4,1:9),
            pos=0,col.axis="red",col.lab="green")
ternaryAxis(side=rep(-1:-3,3),labels=paste(seq(20,80,by=20),"\%"),
            pos=rep(c(0,0.5,1),each=3),col.axis=1:3,col.lab="green")
ternaryAxis(side=rep(1:3,3),at=1:9/10,
            labels=expression(9:1,4:1,7:3,3:2,1:1,2:3,3:7,1:4,1:9),
            pos=rep(c(0,0.5,1),each=3))

plot(acomp(sa.lognormals5),axes=TRUE)
ternaryAxis(side=1:3,pos=0,col.axis="red",col.lab="green")
ternaryAxis(side=1:3,at=1:9/10,
            labels=expression(9:1,4:1,7:3,3:2,1:1,2:3,3:7,1:4,1:9),
            pos=0,col.axis="red",col.lab="green")

}
\keyword{aplot}
