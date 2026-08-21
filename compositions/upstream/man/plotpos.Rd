\name{plot.aplus}
\alias{plot.aplus}
\alias{plot.rplus}
\alias{plot.rmult}
\title{Displaying amounts in scatterplots}
\description{
  This function displays multivariate unclosed amout datasets classes
  "aplus" and "rplus" in a way respecting the choosen geometry
  eventually in log scale. 
}
\usage{
  \method{plot}{aplus}(x,...,labels=colnames(X),cn=colnames(X),
                     aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,
                     center=FALSE,scale=FALSE,pca=FALSE,col.pca=par("col"),
                     add=FALSE,logscale=TRUE,xlim=NULL,ylim=xlim,
                     col=par("col"),plotMissings=TRUE,
                     lenMissingTck=0.05,colMissingTck="red",
                     mp=~simpleMissingSubplot(missingPlotRect,missingInfo,
                                               c("NM","TM",cn)),
                     robust=getOption("robust"))
  \method{plot}{rplus}(x,...,labels=colnames(X),cn=colnames(X),
                     aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,
                     center=FALSE,scale=FALSE,pca=FALSE,col.pca=par("col"),
                     add=FALSE,logscale=FALSE,
                     xlim=NULL,
                     ylim=xlim,col=par("col"),plotMissings=TRUE,
                     lenMissingTck=0.05,colMissingTck="red",
                     mp=~simpleMissingSubplot(missingPlotRect,missingInfo,
                                               c("NM","TM",cn)),
                     robust=getOption("robust"))
  \method{plot}{rmult}(x,...,labels=colnames(X),cn=colnames(X),
                     aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,
                     center=FALSE,scale=FALSE,pca=FALSE,col.pca=par("col"),
                     add=FALSE,logscale=FALSE,col=par("col"),
                     robust=getOption("robust"))
          }
\arguments{
  \item{x}{a dataset with class aplus, rplus or rmult}
  \item{\dots}{further graphical parameters passed (see
    \code{\link[graphics]{par}})}
  \item{add}{a logical indicating whether the information should just
    be added to an existing plot. If FALSE, a new plot is
    created}
  \item{col}{the color to plot the data}
  \item{plotMissings}{ logical indicating that missingness should be
    represented graphically. Componentes with one missing subcomponent
    in the plot are represented by tickmarks at the two
    axis. Cases with two missing components are only
    represented in a special panel drawn according to the \code{mp}
    parameter
    if missings are present. Missings of type BDL (below detection
    limit) are always plotted in nonlogaritmic plots, even if
     \code{plotMissings} is false, but
    in this case this fact is not specially marked.}
  \item{lenMissingTck}{length of the tick-marks (in portion of the
    plotting region) to be plotted for missing
    values. If 0 no tickmarks are plotted. Negative lengths point
    outside of the plot. A length of 1 runs right through the whole plot.}
  \item{colMissingTck}{colors to draw the missing tick-marks. NULL means
  to take the colors specified for the observations.}
  \item{mp}{A formula providing a call to a function plotting
    informations on the missings. The call is evaluted in the
    environment of the panel plotting function and has access (among
    others) to: \code{cn} the names of the components in the current
    plot, \code{x} the dataset of the current plot, \code{missingInfo}
is a table giving the number of 
    observations of the types NM=Non Missing, TM=Totally missing
    (i.e. two components of the subcomposition are missing),
    and the two single component missing possibilities. }

  \item{labels}{the labels for names of the parts}
  \item{cn}{the names of the parts to be used in a single
    panel. Internal use only}
  \item{aspanel}{logical indicating that only a single panel should be
    drawn and not the whole plot. Internal use only}
  \item{id}{a logical. If TRUE one can identify the points like with the
    \code{\link[graphics]{identify}} command}
  \item{idlabs}{A character vector providing the labels to be used with
    the identification, when \code{id=TRUE}}
  \item{idcol}{color of the \code{idlabs} labels}
  \item{center}{a logical indicating whether the data should be
    centered prior to the plot. Centering is done in the chosen
    geometry. See \code{\link{scale}}}
  \item{scale}{a logical indicating whether the data should be
    scaled prior to the plot. Scaling is done in the chosen
    geometry. See \code{\link{scale}}}
  \item{pca}{a logical indicating whether the first principal component
    should be displayed in the plot. Currently, the direction of the
    principal component of the displayed subcomposition is displayed as 
    a line. In a future, the projected principal componenent of the whole 
    dataset should be displayed.}
  \item{col.pca}{the color to draw the principal component.}
  \item{logscale}{logical indicating whether a log scale should be used}
  \item{xlim}{2xncol(x)-matrix giving the xlims for the columns of x}
  \item{ylim}{2xncol(x)-matrix giving the ylims for the columns of x}
  \item{robust}{A robustness description. See
    \link{robustnessInCompositions} for details. The option is used for
    centering, scaling and principle components.}
}
\details{
TO DO: fix pca bug
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{plot.aplus}},
  \code{\link{qqnorm.acomp}},\code{\link{boxplot.acomp}} }

\examples{
data(SimulatedAmounts)
plot(aplus(sa.lognormals))
plot(rplus(sa.lognormals))
plot(aplus(sa.lognormals5))
plot(rplus(sa.lognormals5))
}
\keyword{multivariate}
