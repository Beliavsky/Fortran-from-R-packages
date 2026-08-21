\name{CoDaDendrogram}
\alias{CoDaDendrogram}
\title{Dendrogram representation of acomp or rcomp objects}
\description{
  Function for plotting CoDa-dendrograms of acomp or rcomp objects.
}
\usage{
CoDaDendrogram(X, V = NULL, expr=NULL, mergetree = NULL, signary = NULL, 
    range = c(-4,4), ..., xlim = NULL, ylim = NULL, yaxt = NULL, box.pos = 0,
    box.space = 0.25, col.tree = "black", lty.tree = 1, lwd.tree = 1,
    col.leaf = "black", lty.leaf = 1, lwd.leaf = 1, add = FALSE,border=NULL,
    type = "boxplot")
          }
\arguments{
  \item{X}{ data set to plot (an rcomp or acomp object) }
  \item{V}{ basis to use, described as an ilr matrix }
  \item{expr}{ a formula describing the partition basis, as with \code{\link{balanceBase}} }
  \item{mergetree}{ basis to use, described as a merging tree (as in \code{\link[stats]{hclust}}) }
  \item{signary}{ basis to use, described as a sign matrix (as in the example below)}
  \item{range}{ minimum and maximum value for all coordinates (horizontal axes) }
  \item{...}{further parameters to pass to any function, be it a plotting function or one related to the "type" parameter below; likely to produce lots of warnings}
  \item{xlim}{ minimum and maximum values for the horizontal direction of the plot (related to number of parts) }
  \item{ylim}{ minimum and maximum values for the vertical direction of the plot (related to variance of coordinates) }
  \item{yaxt}{ axis type for the vertical direction of the plot (see \code{\link[graphics]{par}})}
  \item{box.pos}{ if type="boxplot", this is the relative position of the box in the vertical direction: 0 means centered on the axis, -1 aligned below the axis and +1 aligned above the axis }
  \item{box.space}{ if type="boxplot", size of the box in the vertical direction as a portion of the minimal variance of the coordinates }
  \item{col.tree}{ color for the horizontal axes }
  \item{lty.tree}{ line type for the horizontal axes }
  \item{lwd.tree}{ line width for the horizontal axes }
  \item{col.leaf}{ color for the vertical conections between an axis and a part (leaf) }
  \item{lty.leaf}{ line type for the leaves }
  \item{lwd.leaf}{ line width for the leaves }
  \item{add}{ should a new plot be triggered, or is the material to be added to an existing CoDa-dendrogram? }
  \item{border}{the color for drawing the rectangles}
  \item{type}{ what to represent? one of
    "boxplot","density","histogram","lines","nothing" or "points", or an
    univocal abbreviation }
  
}
\details{
  The object \emph{and an isometric basis} are represented in a
  CoDa-dendrogram, as defined by Egozcue and Pawlowsky-Glahn
  (2005). This is a representation of the following elements:
  \describe{
  \item{a}{a hierarchical partition (which can be specified either through an ilrBase matrix (see \code{\link{ilrBase}}), a merging tree structure (see \code{\link[stats]{hclust}}) or a signary matrix (see \code{\link{gsi.merge2signary}}))}
  \item{b}{the sample mean of each coordinate of the ilr basis associated to that partition}
  \item{c}{the sample variance of each coordinate of the ilr basis associated to that partition}
  \item{d}{optionally (potentially!), any graphical representation of each coordinate, as long as this representation is suitable for a univariate data set (box-plot, histogram, dispersion and kernel density are programmed or intended to, but any other may be added with little work).}
}

  Each coordinate is represented in a horizontal axis, which limits correspond to the values given in the parameter \code{range}. The vertical bar going up from each one of these coordinate axes represent the variance of that specific coordinate, and the contact point the coordinate mean. Note that to be able to represent an initial dendrogram, the first call to this function must be given a full data set, as means and variances must be computed. This information is afterwards stored in a global list, to add any sort of new material to all coordinates.
\cr
The default option is \code{type="boxplot"}, which produces a box-plot for each coordinate, customizable using \code{box.pos} and \code{box.space}, as well as typical \code{\link[graphics]{par}} parameters (col, border, lty, lwd, etc.). To obtain only the first three aspects, the function must be called with \code{type="lines"}. As extensions, one might represent a single datum/few data (e.g., a mean or a random subsample of the data set) calling the function with \code{add=TRUE} and \code{type="points"}. Other options (calling functions \code{\link[graphics]{hist}} or \code{\link[stats]{density}}, and admitting their parameters) will be  also soon available.
\cr
Note that the original coda-dendrogram as defined by Egozcue and Pawlowsky-Glahn (2005) works with acomp objects and ilr bases. Functionality is extended to rcomp objects using calls to \code{\link{idt}}.
}

\author{Raimon Tolosana-Delgado, K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{
  \code{\link{ilrBase}},\code{\link{balanceBase}}, \code{\link{rcomp}}, \code{\link{acomp}},
}

\references{
  Egozcue J.J., V. Pawlowsky-Glahn, G. Mateu-Figueras and
  C. Barcel'o-Vidal (2003) Isometric logratio transformations for
  compositional data analysis. \emph{Mathematical Geology}, \bold{35}(3)
  279-300\cr

  Egozcue, J.J. and V. Pawlowsky-Glahn (2005). CoDa-Dendrogram: a new exploratory tool.
  In: Mateu-Figueras, G. and Barcel\'o-Vidal, C. (Eds.)
  \emph{Proceedings of the 2nd International Workshop on Compositional Data Analysis},
  Universitat de Girona, ISBN 84-8458-222-1, \url{https://ima.udg.edu/Activitats/CoDaWork05/}\cr


}
\examples{
# first example: take the data set from the example, select only
# compositional parts
data(Hydrochem)
x = acomp(Hydrochem[,-c(1:5)])
gr = Hydrochem[,4] # river groups (useful afterwards)
# use an ilr basis coming from a clustering of parts
dd = dist(t(clr(x)))
hc1 = hclust(dd,method="ward.D")
plot(hc1)
mergetree=hc1$merge
CoDaDendrogram(X=acomp(x),mergetree=mergetree,col="red",range=c(-8,8),box.space=1)
# add the mean of each river
color=c("green3","red","blue","darkviolet")
aux = sapply(split(x,gr),mean)
aux
CoDaDendrogram(X=acomp(t(aux)),add=TRUE,col=color,type="points",pch=4)

# second example: box-plots by rivers (filled)
CoDaDendrogram(X=acomp(x),mergetree=mergetree,col="black",range=c(-8,8),type="l")
xsplit = split(x,gr)
for(i in 1:4){
 CoDaDendrogram(X=xsplit[[i]],col=color[i],type="box",box.pos=i-2.5,box.space=0.5,add=TRUE)
}

# third example: fewer parts, partition defined by a signary, and empty box-plots
x = acomp(Hydrochem[,c("Na","K","Mg","Ca","Sr","Ba","NH4")])
signary = t(matrix(  c(1,   1,   1,  1,   1,   1,  -1,
                       1,   1,  -1, -1,  -1,  -1,   0,
                       1,  -1,   0,  0,   0,   0,   0,
                       0,   0,  -1,  1,  -1,  -1,   0,
                       0,   0,   1,  0,  -1,   1,   0,
                       0,   0,   1,  0,   0,  -1,   0),ncol=7,nrow=6,byrow=TRUE))

CoDaDendrogram(X=acomp(x),signary=signary,col="black",range=c(-8,8),type="l")
xsplit = split(x,gr)
for(i in 1:4){
  CoDaDendrogram(X=acomp(xsplit[[i]]),border=color[i],
       type="box",box.pos=i-2.5,box.space=1.5,add=TRUE)
  CoDaDendrogram(X=acomp(xsplit[[i]]),col=color[i],
       type="line",add=TRUE)
}
}
\keyword{hplot}
