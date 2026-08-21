\name{CompoundAuxDensPlot}
\Rdversion{2.1}
\alias{CompoundAuxDensPlot}
\alias{dplot.cavgb2}
\title{
Comparison of the compound GB2 and kernel densities by group
}
\description{
Function \code{dplot.cavgb2} produces a plot in which the compound and kernel (Epanechnikov) densities are plotted by group. 
}
\usage{
dplot.cavgb2(group, x, shape1, scale, shape2, shape3, pl0, pl, w=rep(1,length(x)), 
xmax = max(x)*(2/3), ymax=2e-05, decomp="r", choicecol=1:length(levels(group)), 
xlab="") 
}
\arguments{
  \item{group}{numeric; a factor variable giving the group membership of each sampled unit.}
  \item{x}{numeric; can be a vector. The value(s) at which the density is calculated, used for the kernel estimate only. \code{x} is positive. }
  \item{shape1, scale, shape2, shape3}{numeric; positive parameters of the GB2 distribution. On the plot they are denotes as \code{a}, \code{b}, \code{p}, \code{q} and \code{pl0} respectively.}
  \item{pl0}{numeric; a vector of initial proportions defining the number of components and the weight of each component density in the decomposition. Sums to one.}
  \item{pl}{numeric; a vector of fitted proportions (output of \code{\link{pkl.cavgb2}}). Sums to one. If \code{pl} is equal to \code{pl0}, we obtain the GB2 distribution.}
  \item{w}{numeric; weights.}
  \item{xmax}{numeric; scale on the horizontal axis. By default is equal to \eqn{max(x)*(2/3)}.}
  \item{ymax}{numeric; scale on the vertical axis. By default is equal to 2e-05.}
  \item{decomp}{string; specifying if the decomposition of the GB2 is done with respect to the right tail ("r") or the left tail ("l") of the distribution. By default, \code{decomp} = "r" - right tail decomposition.}
  \item{choicecol}{numeric vector of length the number of groups; defines the color with which the density curves will be plotted for each group.}
   \item{xlab}{string; label for \eqn{x}. The default is " ".}
}
\details{The legend is placed interactively.  
}
\value{\code{dplot.cavgb2} plots a graph with two curves - the GB2 density, the compound GB2 per group and the corresponding kernel estimate.  
}
\author{
Monique Graf and Desislava Nedyalkova
}

