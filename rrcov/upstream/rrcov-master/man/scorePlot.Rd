\name{scorePlot-methods}
\docType{methods}
\alias{scorePlot}
\alias{scorePlot-methods}
\alias{scorePlot,ANY-method}
\alias{scorePlot,Pca-method}
\title{
    Score plot for Principal Components (objects of class 'Pca')
}
\description{
  Produces a score plot from an object (derived from)  \code{\link{Pca-class}}.
}
\usage{
    \S4method{scorePlot}{Pca}(x, i=1, j=2, \dots)
}
\arguments{
  \item{x}{an object of class (derived from) \code{"Pca"}.}
  \item{i}{First score coordinate, defaults to \code{i=1}.}
  \item{j}{Second score coordinate, defaults to \code{j=2}.}
  \item{\dots}{optional arguments to be passed to the internal graphical functions.}
}

\section{Side Effects}{
    a plot is produced on the current graphics device.
}
\section{Methods}{
\describe{

    \item{scorePlot}{\code{signature(x = Pca)}: Plot a scatter plot of ith against jth score
    of the Pca object with superimposed tollerance (0.975) ellipse. See also \code{\link{biplot}}, \code{\link{screeplot}}.}
}}

%\references{}
\seealso{
    \code{\link{Pca-class}},
    \code{\link{PcaClassic}},
    \code{\link{PcaRobust-class}}.
}
\examples{
require(graphics)

## PCA of the Hawkins Bradu Kass's Artificial Data
##  using all 4 variables
data(hbk)
pca <- PcaHubert(hbk)
pca

scorePlot(pca)
}
\keyword{multivariate}
\keyword{hplot}
