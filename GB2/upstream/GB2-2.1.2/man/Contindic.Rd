\name{Contindic}
\Rdversion{2.1}
\alias{Contindic}
\alias{contindic.gb2}
\title{
Sensitivity Analysis of Laeken Indicators on GB2 Parameters
}
\usage{
contindic.gb2(resol, shape1, shape21, shape22, shape31, shape32, fn, title, table=FALSE)
}
\arguments{
  \item{resol}{numeric; number of grid points horizontally and vertically.}
  \item{shape1}{numeric; positive parameter, first shape parameter of the GB2 distribution.} 
  \item{shape21, shape22, shape31, shape32}{numeric; limits on the positive parameters of the Beta distribution.}
  \item{fn}{string; the name of the function to be used for the calculation of the values to be plotted.}
  \item{title}{string; title of the plot.}
  \item{table}{boolean; if \code{TRUE}, a table containing the values of the function \code{fn} at the different grid points is printed.}
}
\description{Produces a contour plot of an indicator for a given \code{shape1}.
}
\details{An indicator is defined as a function of three parameters. The shape parameter, \code{shape1}, is held fixed. The shape parameters \code{shape2} and \code{shape3} vary between \code{shape21} and \code{shape22}, and \code{shape31} and \code{shape32}, respectively.
}
\value{
A contour plot of a given indicator for a fixed value of the shape parameter \code{shape1}.
}
\seealso{\code{\link{contour}} (package \code{graphics}) for more details on contour plots.
}
\author{
Monique Graf
}
\examples{
par(mfrow=c(2,2))
shape21 <- 0.3
shape31 <- 0.36
shape22 <- 1.5
shape32 <- 1.5
shape11 <- 2.7
shape12 <- 9.2
resol <- 11
rangea <- round(seq(shape11, shape12 ,length.out=4),digits=1)
arpr <- function(shape1, shape2, shape3) 100*arpr.gb2(0.6, shape1, shape2, shape3)
fonc <- "arpr"
for (shape1 in rangea){
contindic.gb2(resol, shape1, shape21, shape22, shape31, shape32, arpr, "At-risk-of-poverty rate", 
table=TRUE)
}
}

