\name{ilt}
\alias{ilt}
\alias{iltInv}
\title{Isometric log transform}
\description{
  Compute the isometric log  transform of a vector (dataset) of
  amounts and its inverse.
}
\usage{
          ilt( x ,...)
          iltInv( z ,... )
          }
\arguments{
  \item{x}{a vector or data matrix of amounts}
  \item{z}{the ilt-transform of a vector or data matrix of
    ilt-transforms of amounts}
  \item{\dots}{generic arguments, not used.}
}
\value{
  \code{ilt} gives the isometric log transform, i.e. simply the log of
  the argument, whereas
  \code{iltInv} gives amounts with the given ilt, i.e. simply the exp
  of the argument.
  }
\details{
  The ilt-transform maps  D amounts (considered in log geometry)
  isometrically to a D dimensional euclidean vector. The \code{ilt} is
  part of the \code{\link{aplus}} framework.
  \cr 
  The data can then be analysed in this transformation by all classical
  multivariate analysis tools. The interpretation of the results is easy
  since the relation to the original
  variables is preserved. \cr
    
  The isometric log transform is given by
    \deqn{ ilt(x)_i := \ln x_i }

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{ilr}}, \code{\link{iit}}, \code{\link{aplus}}
}
\references{
    van den Boogaart, K.G. and R. Tolosana-Delgado (2008) "compositions": a unified 
  R package to analyze Compositional Data, \emph{Computers &
    Geosciences}, 34 (4), pages 320-338, \doi{10.1016/j.cageo.2006.11.017}.
}
\examples{
(tmp <- ilt(c(1,2,3)))
iltInv(tmp)
iltInv(tmp) - c(1,2,3) # 0
data(Hydrochem)
cdata <- Hydrochem[,6:19]
pairs(ilt(cdata)) 
}
\keyword{multivariate}
