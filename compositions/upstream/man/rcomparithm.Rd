\name{rcomparithm}
\alias{+.rcomp}
\alias{-.rcomp}
\alias{*.rcomp}
\alias{/.rcomp}
\alias{convex.rcomp}
\title{Arithmetic operations for compositions in a real geometry}
\description{
  The real compositions form a manifold of the real vector space. The
  induced operations +,-,*,/ give results valued in the real 
  vector space, but possibly outside the simplex. 
}
\usage{
convex.rcomp(x,y,alpha=0.5)
## Methods for class "rcomp"
##   x+y
##   x-y
##   -x
##   x*r
##   r*x
##   x/r
          }
\arguments{
  \item{x}{an rcomp composition or dataset of compositions }
  \item{y}{an rcomp composition or dataset of compositions }
  \item{r}{a numeric vector of size 1 or nrow(x)}
  \item{alpha}{a numeric vector of size 1 or nrow(x) with values between
  0 and 1}
}
\value{
  \code{rmult}-objects containing the given operations on the simplex
   as subset of the \eqn{R^D}. Only the convex combination
  \code{convex.rcomp} results in an \code{rcomp}-object again, since
  only this operation is closed.
}
\details{
  The functions behave quite like \code{\link{+.rmult}}.\cr
  The convex combination is defined as: \code{x*alpha + (1-alpha)*y}
}
\note{For \code{*} the arguments x and y can be exchanged.}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{+.rmult}}, \code{\link{+.acomp}},\code{\link{cpt}}, \code{\link{rcomp}}, \code{\link{rmult}}
}

\examples{
rcomp(1:5)* -1 + rcomp(1:5)
data(SimulatedAmounts)
cdata <- rcomp(sa.lognormals)
plot( tmp <- (cdata-mean(cdata))/msd(cdata) )
class(tmp)
mean(tmp)
msd(tmp)
var(tmp)
plot(convex.rcomp(rcomp(c(1,1,1)),sa.lognormals,0.1))
}
\keyword{multivariate}
