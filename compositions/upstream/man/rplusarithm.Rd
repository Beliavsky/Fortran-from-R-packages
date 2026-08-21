\name{rplusarithm}
\alias{+.rplus}
\alias{-.rplus}
\alias{*.rplus}
\alias{/.rplus}
\alias{mul.rplus}
\title{vectorial arithmetic for data sets with rplus class}
\description{
  The positive quadrant forms a manifold of the real vector space. The
  induced operations +,-,*,/ give results valued in this real vector space 
  (not necessarily inside the manifold).
}
\usage{
mul.rplus(x,r)
## Methods for class rplus
##   x+y
##  x-y
##   -x
##   x*r
##   r*x
##   x/r
          }
\arguments{
  \item{x}{an rplus composition or dataset of compositions }
  \item{y}{an rplus composition or dataset of compositions }
  \item{r}{a numeric vector of size 1 or nrow(x)}
}
\value{
  \code{rmult}-objects containing the given operations on the rcomp
  manifold as subset of the \eqn{R^D}. Only the addition and
  multiplication with positive numbers are internal
  operation and results in an \code{rplus}-object again.
}
\details{
  The functions behave quite like \code{\link{+.rmult}}.\cr
}
\note{For \code{*} the arguments x and y can be exchanged.}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{+.rmult}}, \code{\link{+.acomp}},\code{\link{cpt}}, \code{\link{rcomp}}, \code{\link{rmult}}
}

\examples{
rplus(1:5)* -1 + rplus(1:5)
data(SimulatedAmounts)
cdata <- rplus(sa.lognormals)
plot( tmp <- (cdata-mean(cdata))/msd(cdata) )
class(tmp)
mean(tmp)
msd(tmp)
var(tmp)
}
\keyword{multivariate}
