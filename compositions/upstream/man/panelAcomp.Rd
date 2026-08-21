\name{pairs}
\alias{pairs.acomp}
\alias{pairs.rcomp}
\alias{vp.lrdensityplot}
\alias{vp.diffdensityplot}
\alias{vp.lrdensityplot}
\alias{vp.lrboxplot}
\alias{vp.kde2dplot}
\title{Pairs plot method for compositions}
\description{
  Pairs plot function for compositions, allowing flexible representations.
}
\usage{
 \method{pairs}{acomp}(x, labels, panel = vp.lrdensityplot, ...,
                         horInd = 1:ncol(x), verInd = 1:ncol(x),
                         lower.panel = panel, upper.panel = panel,
                         diag.panel = NULL, text.panel = textPanel,
                         label.pos = 0.5 + has.diag/3, line.main = 3,
                         cex.labels = NULL, font.labels = 1,
                         row1attop = TRUE, gap = 1, log = "")
 \method{pairs}{rcomp}(x, labels, panel = vp.diffdensityplot, ...,
                         horInd = 1:ncol(x), verInd = 1:ncol(x),
                         lower.panel = panel, upper.panel = panel,
                         diag.panel = NULL, text.panel = textPanel,
                         label.pos = 0.5 + has.diag/3, line.main = 3,
                         cex.labels = NULL, font.labels = 1,
                         row1attop = TRUE, gap = 1, log = "")
  vp.lrdensityplot(x, y, col=2,..., alpha = NULL)
  vp.diffdensityplot(x, y, col=2,..., alpha = NULL)
  vp.lrboxplot(x, y, ...)
  vp.kde2dplot(x, y, grid=TRUE, legpos="bottomright", colpalette=heat.colors,...)
}
\arguments{
  \item{x}{a dataset of a compositional class; 
  or for the panel functions, a vector of row components}
  \item{y}{for the panel functions, a vector of column components}
  \item{\dots}{further graphical parameters passed (see
    \code{\link[graphics]{par}})}
  \item{labels}{the names of the parts}
  \item{panel}{common panel function to use for all off-diagonal plots}
  \item{horInd}{indices of columns of x to plot on the horizontal axis, defaults to all columns} 
  \item{verInd}{indices of columns of x to plot on the vertical axis, defaults to all columns}
  \item{lower.panel}{panel function for the lower triangle of plots, defaults to the common panel}
  \item{upper.panel}{panel function for the uppper triangle of plots, defaults to the common panel}
  \item{diag.panel}{panel function for the diagonal of plots, defaults to text.panel}
  \item{text.panel}{panel function to write labels on the diagonal panels}
  \item{label.pos}{y position of labels in the text panel} 
  \item{line.main}{if main is specified, line.main gives the line argument to mtext() which draws the title. You may want to specify oma when changing line.main}
  \item{cex.labels}{graphics parameters for the text panel}
  \item{font.labels}{graphics parameters for the text panel}
  \item{row1attop}{logical. Should the layout be matrix-like with row 1 at the top, or graph-like with row 1 at the bottom?} 
  \item{gap}{distance between subplots, in margin lines} 
  \item{log}{a character string indicating if logarithmic axes are to be used: see plot.default. Should not be used and left to the panel function to handle}
  \item{col}{color for density and histogram components of the panel vp.*density} 
  \item{alpha}{alpha level for marking normality in the panels vp.*density; default to no mark} 
  \item{grid}{should a unit-grid be added to each panel?}
  \item{legpos}{where should the legend be placed? to be given as \code{x} in \code{graphics::legend}}
  \item{colpalette}{which color palette is desired for the 2d density levels?}
}
\details{
  The data is displayed in a matrix of plots, after the indications of a panel function. 
  This is a simple implementation of \code{\link{pairs}} compositional methods, the real 
  functionality is controlled by the panel functions. 
  
  The three panel functions included here can be used for generating either boxplots
  or histograms plus kernel density plots of all pairwise logratios (in acomp) or 
  differences (in rcomp) of the components. In the cas of histograms, these 
  can be colored or left black-and-white depending on the adjustment to
  normality, controlled by a \code{\link{shapiro.test}} and the alpha-level given.
  These panel functions serve also as examples of how to generate user defined panels. 
}

\author{Raimon Tolosana-Delgado}

\seealso{
  \code{\link{vp.boxplot}}, \code{\link{vp.logboxplot}}
}

\examples{
data(SimulatedAmounts)
pairs(acomp(sa.lognormals))
pairs(rcomp(sa.lognormals))
}

\keyword{hplot}
