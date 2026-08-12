\name{fminviz,flineviz}
\alias{fminviz}
\alias{flineviz}
\title{
  Visualize Function Minimum
}
\description{
  Visualizes multivariate functions around a point or along a line between 
  two points in \code{R^n}.
}
\usage{
fminviz(fn, x0, nlines = 2*length(x0),
        npoints = 51, scaled = 1.0)

flineviz(fn, x1, x2, npoints = 51, scaled = 0.1)
}
\arguments{
  \item{fn}{multivariate function to be visualized.}
  \item{x0,x1,x2}{points in n-dimensional space.}
  \item{nlines}{number of lines to plot.}
  \item{npoints}{number of points used to plot a line.}
  \item{scaled}{scale factor to extend the line(s).}
}
\details{
  \code{fminviz} vizualizes the behavior of a multivariate function \code{fn} 
  around a point \code{x0}. It randomly selects \code{nlines} lines through 
  \code{x0} in \code{R^n} and draws the curves of the function along these 
  lines in one graph.

  Curves that have at least one point below \code{fn(x0)} are drawn in red, 
  all others in blue. The scale on the x-axis is the Euclidean distance in 
  \code{R^n}. The \code{scale} factor can change it.

  \code{flineviz} vizualizes the behavior of a multivariate function \code{fn} 
  along the straight line between the points \code{x1} and \code{x2}. Points
  \code{x1} and \code{x2} are also plotted.
}
\value{
  Plots a line graph and returns \code{NULL} (invisibly).
}
\examples{
\dontrun{
  f1 <- function(x) x[1]^2 - x[2]^2
  fminviz(f1, c(0, 0), nlines = 10)

  f2 <- function(x) (1 - x[1])^2 + 100*(x[2] - x[1]^2)^2
  flineviz(f2, c(0, 0), c(1, 1))
  }
}

\keyword{ visualize }
