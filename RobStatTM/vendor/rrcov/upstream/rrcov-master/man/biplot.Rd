\name{biplot-methods}
\docType{methods}
\alias{biplot}
\alias{biplot-methods}
\alias{biplot,ANY-method}
\alias{biplot,Pca-method}
\title{
    Biplot for Principal Components (objects of class 'Pca')
}
\description{
  Produces a biplot from an object (derived from)  \code{\link{Pca-class}}.
}
\usage{
    \S4method{biplot}{Pca}(x, choices=1L:2L, scale=1, \dots)
}
\arguments{
  \item{x}{an object of class (derived from) \code{"Pca"}.}
  \item{choices}{length 2 vector specifying the components to plot. Only the default is a biplot in the strict sense.}
  \item{scale}{
    The variables are scaled by \code{lambda ^ scale} and the
    observations are scaled by \code{lambda ^ (1-scale)} where
    \code{lambda} are the singular values as computed by
    the Principal Components function. Normally \code{0 <= scale <= 1}, and a warning
    will be issued if the specified \code{scale} is outside this range.
  }
  \item{\dots}{optional arguments to be passed to the internal graphical functions.}
}

\section{Side Effects}{
    a plot is produced on the current graphics device.
}
\section{Methods}{
\describe{

    \item{biplot}{\code{signature(x = Pca)}: Plot a biplot, i.e. represent both
        the observations and variables of a matrix of multivariate data on the same
        plot. See also \code{\link{biplot.princomp}}.}
}}

\references{
    Gabriel, K. R. (1971).
    The biplot graphical display of matrices with applications to
    principal component analysis.
    \emph{Biometrika}, \bold{58}, 453--467.
}
\seealso{
    \code{\link{Pca-class}},
    \code{\link{PcaClassic}},
    \code{\link{PcaRobust-class}}.
}
\examples{
require(graphics)
biplot(PcaClassic(USArrests, k=2))
}
\keyword{multivariate}
\keyword{hplot}
