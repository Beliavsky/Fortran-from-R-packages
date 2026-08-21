\name{plot.acomp}
\alias{plot.acomp}
\alias{plot.rcomp}
\alias{plot.ccomp}
\title{Ternary diagrams}
\description{
  Displaying compositions in ternary diagrams
}
\usage{
 \method{plot}{acomp}(x,...,labels=names(x),
          aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,center=FALSE,
          scale=FALSE,pca=FALSE,col.pca=par("col"),margin="acomp",
          add=FALSE,triangle=!add,col=par("col"),axes=FALSE,
          plotMissings=TRUE,
          lenMissingTck=0.05,colMissingTck="red",
          mp=~simpleMissingSubplot(c(0,1,0.95,1),
                                  missingInfo,c("NM","TM",cn)),
          robust=getOption("robust"))
 \method{plot}{rcomp}(x,...,labels=names(x),
          aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,center=FALSE,
          scale=FALSE,pca=FALSE,col.pca=par("col"),margin="rcomp",
          add=FALSE,triangle=!add,col=par("col"),axes=FALSE
          ,plotMissings=TRUE,
          lenMissingTck=0.05,colMissingTck="red",
          mp=~simpleMissingSubplot(c(0,1,0.95,1),
                                   missingInfo,c("NM","TM",cn)),
          robust=getOption("robust"))
 \method{plot}{ccomp}(x,...)
          }
\arguments{
  \item{x}{a dataset of a compositional class}
  \item{\dots}{further graphical parameters passed (see
    \code{\link[graphics]{par}})}
  \item{margin}{the type of marginalisation to be computed, when
    displaying the individual panels. Possible values are: \code{"acomp"},
    \code{"rcomp"} and any of the variable names/column numbers in the
    composition. If one of the columns is selected each panel displays a
    subcomposition given by the row part, the column part and
    the given part. If one of the classes is given the corresponding
    margin \code{\link{acompmargin}}  or  \code{\link{rcompmargin}} is
    used. }
  \item{add}{a logical indicating whether the information should just
    be added to an existing plot. If FALSE a new plot is created}
  \item{triangle}{a logical indicating whether the triangle should be
    drawn}
  \item{col}{the color to plot the data}
  \item{labels}{the names of the parts}
%  \item{cn}{the names of the parts to be used in a single
%    panel. Internal use only.}
  \item{aspanel}{logical indicating that only a single panel should be
    drawn and not the whole plot. Internal use only}
  \item{id}{logical, if TRUE one can identify the points like with the
    \code{\link[graphics]{identify}} command. }
  \item{idlabs}{a character vector providing the labels to be used with
    the identification, when \code{id=TRUE}}
  \item{idcol}{color of the \code{idlabs} labels}
  \item{center}{a logical indicating whether a the data should be
    centered prior to the plot. Centering is done in the choosen
    geometry. See \code{\link{scale}}}
  \item{scale}{a logical indicating whether a the data should be
    scaled prior to the plot. Scaling is done in the choosen
    geometry. See \code{\link{scale}}}
  \item{pca}{a logical indicating whether the first principal component
    should be displayed in the plot. Currently, the direction of the
    principal component of the displayed subcomposition is displayed as 
    a line. In a future, the projected principal componenent of the whole 
    dataset should be displayed.}
  \item{col.pca}{The color to draw the principal component.}
  \item{axes}{Either a logical wether to plot the axes, or numerical
    enumerating the axes sides to be used e.g. 1 for only plotting the
    lower axes, or a list of parameters to ternaryAxis. }
  \item{plotMissings}{ logical indicating that missingness should be
    represented graphically. Componentes with one missing subcomponent
    in the plot are represented by tickmarks at the three
    axis. Components with two or three missing components are only
    represented in a special panel drawn according to the mp parameter
    if missings are present. Missings of type BDL (below detection
    limit) are always plotted, even if \code{plotMissings} is false, but
    in this case this fact is not specially marked. In rcomp geometry an
    actuall 0 in the data is never treated as missing.}
  \item{lenMissingTck}{length of the tick-marks to be plotted for missing
    values. If 0 no tickmarks are plotted. Negative lengths point
    outside. length 1 draws right through to the opposit
    corner. Missing ticks in acomp geometry are inclined showing the
    line of possible values in acomp geometry. Missingticks in
    rcomp-geometry are vertical to the axis representing the fact that
    only the other component is unkown. That these lines can leave the plot
    is one of the  odd consequences of rcomp geometry.}
  \item{colMissingTck}{colors to draw the missing tick-marks. NULL means
  to take the colors specified for the observations.}
  \item{mp}{A formula providing a call to a function plotting
    informations on the missings. The call is evaluted in the
    environment of the panel plotting function and has access (among
    others) to: \code{cn} the names of the components in the current
    plot, \code{x} the dataset of the current plot, \code{y} the
    transformed dataset, (c60,s60) coordinates of the upper vertex of
    the triangle. \code{missingInfo} is a table giving the number of
    observations of the types NM=Non Missing, TM=Totally missing
    (i.e. at least two components of the subcomposition are missing),
    and the three single component missing possibilities for the three
    components. }
  \item{robust}{A robustness description. See
    \link{robustnessInCompositions} for details. The option is used for
    centering, scaling and principle components. }
}
\details{
  The data is displayed in ternary diagrams. Thus, it does not work for
  two-part compositions. Compositions of three parts are displayed
  in a single ternary diagram. For compositions of more than three
  components, the data is arranged in a scatterplot matrix through the
  command \code{\link[graphics]{pairs}}.
  \cr
  In this case, the third component in each of the panels is chosen 
  according to setting of \code{margin=}. Possible values of \code{margin=} are:
  \code{"acomp"}, \code{"rcomp"} and any of the variable names/column numbers in the
  composition. If one of the columns is selected each panel displays a
  subcomposition given by the row part, the column part and
  the given part. If one of the classes is given the corresponding
  margin \code{\link{acompmargin}}  or  \code{\link{rcompmargin}} is
  used.
  \cr
  Ternary diagrams can be read in multiple ways. Each corner of the
  triangle corresponds to an extreme composition containing only the part
  displayed in that corner. Points on the edges correspond to
  compositions containing only the parts in the adjacent corners. The
  relative amounts are displayed by the distance to the opposite
  corner (so-called barycentric coordinates). The individual portions 
  of any point can be infered by drawing a line through the investigated point, 
  and parallel to the edge opposite to the corner of the part of interest. 
  The portion of this part is constant along the line. Thus we can read it 
  on the sides of the ternary diagram, where the line crosses its borders.
  Note that these \code{\link{isoPortionLines}} remain straight under an 
  arbitrary perturbation. 
  \cr
  ccomp ternary diagrams are always jittered to avoid overplotting.
} \references{
    Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
      Data} Monographs on Statistics and Applied Probability. Chapman &
    Hall Ltd., London (UK). 416p.\cr

  Aitchison, J, C. Barcel'o-Vidal, J.J. Egozcue, V. Pawlowsky-Glahn
  (2002) A consise guide to the algebraic geometric structure of the
  simplex, the sample space for compositional data analysis, \emph{Terra
    Nostra}, Schriften der Alfred Wegener-Stiftung, 03/2003\cr

  Billheimer, D., P. Guttorp, W.F. and Fagan (2001) Statistical interpretation of species composition,
  \emph{Journal of the American Statistical Association}, \bold{96} (456), 1205-1214\cr

  Pawlowsky-Glahn, V. and J.J. Egozcue (2001) Geometric approach to
  statistical analysis on the simplex. \emph{SERRA} \bold{15}(5), 384-398\cr
  
  \url{https://ima.udg.edu/Activitats/CoDaWork03/}\cr
  
  \url{https://ima.udg.edu/Activitats/CoDaWork05/}

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}

\seealso{
  \code{\link{plot.aplus}}, \code{\link{plot3D}} (for 3D plot),
   \code{\link{kingTetrahedron}} (for 3D-plot model export),
  \code{\link{qqnorm.acomp}},\code{\link{boxplot.acomp}} }

\examples{
data(SimulatedAmounts)
plot(acomp(sa.lognormals))
plot(acomp(sa.lognormals),axes=TRUE)
plot(rcomp(sa.lognormals))
plot(rcomp(sa.lognormals5))

plot(acomp(sa.lognormals5),pca=TRUE,col.pca="red")
plot(rcomp(sa.lognormals5),pca=TRUE,col.pca="red",axes=TRUE)
}
\keyword{hplot}
