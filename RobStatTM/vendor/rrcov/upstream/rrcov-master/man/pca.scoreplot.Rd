\name{pca.scoreplot}
\alias{pca.scoreplot}
\title{
    Score plot for Principal Components (objects of class 'Pca')
}
\description{
  Produces a score plot from an object (derived from)  \code{\link{Pca-class}}.
}
\usage{
    pca.scoreplot(obj, i=1, j=2, main, id.n, \dots)
}

\arguments{
  \item{obj}{an object of class (derived from) \code{"Pca"}.}
  \item{i}{First score coordinate, defaults to \code{i=1}.}
  \item{j}{Second score coordinate, defaults to \code{j=2}.}
  \item{main}{The main title of the plot.}
  \item{id.n}{Number of observations to identify by a label. If missing and the
    total number of observations is less or equal to 10, all observations will
    be labelled.}
  \item{\dots}{Optional arguments to be passed to the internal graphical functions.}
}
%\details{}
%\value{}
%\references{}
%\note{}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
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
pca.scoreplot(pca)
}
\keyword{robust}
\keyword{multivariate}
