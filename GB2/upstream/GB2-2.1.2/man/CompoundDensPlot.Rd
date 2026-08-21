\name{CompoundDensPlot}
\Rdversion{2.1}
\alias{CompoundDensPlot}
\alias{dplot.cgb2}
\title{
Comparison of the GB2, compound GB2 and kernel densities
}
\description{
Function \code{dplot.cgb2} produces a plot in which the three densities are plotted. 
}
\usage{
dplot.cgb2(x,shape1, scale, shape2, shape3, pl0, pl, w=rep(1,length(x)), decomp="r", 
xmax = max(x)*(2/3), choicecol=1:3, kernel="epanechnikov", adjust=1, title=NULL,
ylim=NULL) 
}
\arguments{
  \item{x}{numeric; can be a vector. The value(s) at which the density is calculated, used for the kernel estimate only. \code{x} is positive. }
  \item{shape1, scale ,shape2, shape3}{numeric; positive parameters of the GB2 distribution. On the plot they are denoted as \code{a}, \code{b}, \code{p}, \code{q} and \code{pl0} respectively.}
  \item{pl0}{numeric; a vector of initial proportions defining the number of components and the weight of each component density in the decomposition. Sums to one.}
  \item{pl}{numeric; a vector of mixture probabilities (output of \code{\link{ml.cgb2}}). Sums to one. If \code{pl} is equal to \code{pl0}, we obtain the GB2 distribution.}
  \item{w}{numeric; weights.}
  \item{decomp}{string; specifying if the decomposition of the GB2 is done with respect to the right tail ("r") or the left tail ("l") of the distribution. By default, \code{decomp} = "r" - right tail decomposition.}
  \item{xmax}{numeric; maximum \code{x} value to be plotted.}
  \item{choicecol}{numeric vector of length 3; defines the color with which the density curves will be plotted.}
  \item{adjust}{numeric; graphical parameter of the generic function \code{\link{density}}.}
  \item{title}{string; title of the plot. By default is equall to NULL (no title).}
  \item{ylim}{string; scaling of parameters. By default is equall to NULL (automatic scaling).}
  \item{kernel}{string; the kernel used for the kernel density estimate. The default value is  "Epanechnikov" (see \code{\link{plot.density}}).}
}
\details{ The legend is placed interactively.  
}
\value{\code{dplot.cgb2} plots a graph with three curves - the GB2 density, the compound GB2 density and the corresponding kernel estimate  
}
\author{
Monique Graf and Desislava Nedyalkova
}

