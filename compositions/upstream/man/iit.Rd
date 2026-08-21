\name{iit}
\alias{iit}
\alias{iitInv}
\title{Isometric identity transform}
\description{
  Compute the isometric identity  transform of a vector (dataset) of
  amounts and its inverse.
}
\usage{
          iit( x ,...)
          iitInv( z ,... )
          }
\arguments{
  \item{x}{a vector or data matrix of amounts}
  \item{z}{the iit-transform of a vector or  data.matrix of
    iit-transforms of amounts}
  \item{\dots}{generic arguments, to pass to other functions.}
}
\value{
  \code{ilt} gives the isometric identity transform, i.e. simply the
  input stripped of the "rplus" class attribute,
  \code{iptInv} gives amounts with class "rplus" with the given iit,
  i.e. simply the argument checked to be a valid "rplus" object, and 
  with this class attribute.
  }
\details{
  The iit-transform maps D amounts (considered in a real geometry)
  isometrically to a D dimensonal euclidian vector. The \code{iit} is
  part of the \code{\link{rplus}} framework. Despite its trivial
  operation, it is present to achieve maximal analogy between the
  \code{\link{aplus}} and the \code{\link{rplus}} framework.
  \cr 
  The data can then be analysed in this transformated space by all classical
  multivariate analysis tools. The interpretation of the results is easy
  since the relation to the original
  variables is preserved. However results may be inconsistent, since the
  multivariate analysis tools disregard the positivity condition and the
  inner laws of amounts.\cr
    
  The isometric identity transform is a simple identity given by
    \deqn{ iit(x)_i :=  x_i }

}
\note{\code{iit} can be used to unclass amounts.}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{ilt}}, \code{\link{ilr},} \code{\link{rplus}}
}
\references{
      van den Boogaart, K.G. and R. Tolosana-Delgado (2008) "compositions": a unified 
  R package to analyze Compositional Data, \emph{Computers &
    Geosciences}, 34 (4), pages 320-338, \doi{10.1016/j.cageo.2006.11.017}.

}
\examples{
(tmp <- iit(c(1,2,3)))
iitInv(tmp)
iitInv(tmp) - c(1,2,3) # 0
data(Hydrochem)
cdata <- Hydrochem[,6:19]
pairs(iit(cdata)) 
}
\keyword{multivariate}
