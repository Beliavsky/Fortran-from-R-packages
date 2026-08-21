\name{boxplot}
\alias{boxplot.acomp}
\alias{boxplot.rcomp}
\alias{boxplot.rplus}
\alias{boxplot.aplus}
\alias{vp.boxplot}
\alias{vp.logboxplot}
\title{Displaying compositions and amounts with box-plots}
\description{
  For the different interpretations of amounts or compositional data, a different type of
  boxplot is feasible. Thus different boxplots are drawn.
}
\usage{
\method{boxplot}{acomp}(x,fak=NULL,\dots,
                         xlim=NULL,ylim=NULL,log=TRUE,
                         panel=vp.logboxplot,dots=!boxes,boxes=TRUE,
                          notch=FALSE,
                          plotMissings=TRUE,
                          mp=~simpleMissingSubplot(missingPlotRect,
                                                missingInfo,c("NM","TM",cn))
                          )
\method{boxplot}{rcomp}(x,fak=NULL,\dots,
                         xlim=NULL,ylim=NULL,log=FALSE,
                         panel=vp.boxplot,dots=!boxes,boxes=TRUE,
                          notch=FALSE,
                          plotMissings=TRUE,
                          mp=~simpleMissingSubplot(missingPlotRect,
                                                missingInfo,c("NM","TM",cn)))
\method{boxplot}{aplus}(x,fak=NULL,\dots,log=TRUE,
                          plotMissings=TRUE,
                          mp=~simpleMissingSubplot(missingPlotRect,
                                                   missingInfo,
                                                   names(missingInfo)))
\method{boxplot}{rplus}(x,fak=NULL,\dots,ylim=NULL,log=FALSE,
                          plotMissings=TRUE,
                          mp=~simpleMissingSubplot(missingPlotRect,
                                                   missingInfo,
                                                   names(missingInfo)))
vp.boxplot(x,y,...,dots=FALSE,boxes=TRUE,xlim=NULL,ylim=NULL,log=FALSE,
                          notch=FALSE,plotMissings=TRUE,
                          mp=~simpleMissingSubplot(missingPlotRect,
                                                   missingInfo,c("NM","TM",cn)),
                          missingness=attr(y,"missingness") ) 
vp.logboxplot(x,y,...,dots=FALSE,boxes=TRUE,xlim,ylim,log=TRUE,notch=FALSE,
                          plotMissings=TRUE, 
                          mp=~simpleMissingSubplot(missingPlotRect,
                                                   missingInfo,c("NM","TM",cn)),
                          missingness=attr(y,"missingness")) 
}
\arguments{
  \item{x}{a data set}
  \item{fak}{a factor to split the data set, not yet implemented in aplus
  and rplus}
  \item{xlim}{x-limits of the plot.}
  \item{ylim}{y-limits of the plot.}
  \item{log}{logical indicating whether ploting should be done on log
    scale}
  \item{panel}{the panel function to be used or a list of multiple panel
    functions}
  \item{\dots}{further graphical parameters}
  \item{dots}{a logical indicating whether the points should be drawn}  
  \item{boxes}{a logical indicating whether the boxes should be drawn}
  \item{y}{used by pairs}
  \item{notch}{logical, should the boxes be notched?}
  \item{plotMissings}{Logical indicating that missings should be
displayed.}
\item{mp}{A formula providing a function call, which will be evaluated
within each panel with missings to plot the missingness situation. The
call can use the variables \code{missingPlotRect}, which provides a
rectangle to plot the information to in a par("usr") like
specification. In the r\code{X} is the current data }
\item{missingness}{The missingness information as a result from
\code{\link{missingType}} of the full data information the panels could base there
missing plots on.}
}
\details{
  \code{boxplot.aplus} and \code{boxplot.rplus} are wrappers of \code{\link{bxp}}, 
  which just take into account the possible logarithmic scale of the data. \cr

  \code{boxplot.acomp} and \code{boxplot.rcomp} generate a matrix of box-plots, where 
  each cell represents the difference between the row and column variables. Such 
  difference is respectively computed as a log-ratio and a rest.\cr

  \code{vp.boxplot} and \code{vp.logboxplot} are only used as panel functions. 
  They should not be directly called.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{plot.acomp}}, \code{\link{qqnorm.acomp}} 
}

\examples{
data(SimulatedAmounts)
boxplot(acomp(sa.lognormals))
boxplot(rcomp(sa.lognormals))
boxplot(aplus(sa.lognormals))
boxplot(rplus(sa.lognormals))
# And now with missing!!!
boxplot(acomp(sa.tnormals))

}
\keyword{multivariate}
