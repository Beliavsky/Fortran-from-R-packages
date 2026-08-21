\name{MahalanobisDist}
\alias{MahalanobisDist}
\alias{MahalanobisDist.rmult}
\alias{MahalanobisDist.acomp}
\alias{MahalanobisDist.rcomp}
\alias{MahalanobisDist.rplus}
\alias{MahalanobisDist.aplus}
\title{Compute Mahalanobis distances based von robust Estimations}
\description{
  MahalanobisDist computes the Mahalanobis distances to the center or to
  other observations. 
}
\usage{
MahalanobisDist(x,center=NULL,cov=NULL,inverted=FALSE,...)
\method{MahalanobisDist}{rmult}(x,center=NULL,cov=NULL,inverted=FALSE,...,
           goodOnly=NULL,pairwise=FALSE,pow=1,
robust=FALSE,giveGeometry=FALSE)
\method{MahalanobisDist}{acomp}(x,center=NULL,cov=NULL,inverted=FALSE,...,
           goodOnly=NULL, pairwise=FALSE,pow=1,robust=FALSE,giveGeometry=FALSE)
          }
\arguments{
  \item{x}{the dataset}
  \item{robust}{logical or a robust method description (see
    \code{\link{robustnessInCompositions}}) specifiying how the center and covariance
    matrix are estimated,if not given.}
  \item{\dots}{Further arguments to \code{\link{solve}}. }
  \item{center}{An estimated for the center (mean) of the dataset. If
    center is NULL it will be estimated based using the given robust
    option.}
  \item{cov}{An estimated for the spread (covariance matrix) of the
dataset. If
    cov is NULL it will be estimated based using the given robust
    option.}
  \item{inverted}{TRUE if the inverse of the covariance matrix is
    given.}
  \item{goodOnly}{An vector of indices to the columns of x that should
    be used for estimation of center and spread.}
  \item{pairwise}{If FALSE the distances to the center are returned as a
  vector. If TRUE the distances between the cases are returned as a
  distance matrix.}
   \item{pow}{The power of the Mahalanobis distance to be used. 1
     correponds to the square root of the squared distance in
     transformed space, like it is defined in most books. The choice 2
     corresponds to what is implemented in many software package
     including the \code{\link{mahalanobis}} function in R.}
  \item{giveGeometry}{If true an atrributes \code{"center"} and
    \code{"cov"} given the center and the idt-variance used for the
    calculations. 
  }
}
\value{
  Either a vector of Mahalanobis distances to the center, or a distance
  matrix (like from \code{\link{dist}}) giving the pairwise Mahalanobis
  distances of the data.
}
\details{
  The Mahalanobis distance is the distance in a linearly transformed
  space, where the linear transformation is selected in such a way,that
  the variance is the unit matrix. Thus the distances are given in
  multiples of standard deviation. 
}
\note{
  Unlike the \code{\link{mahalanobis}} function this function does not
  be default compute the square of the mahalanobis distance. The pow
  option is provided if the square is needed.
  \cr
  The package \pkg{robustbase} is required for using the
  robust estimations.}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{dist}}, \code{\link{OutlierClassifier1}}
}

\examples{
data(SimulatedAmounts)
data5 <- acomp(sa.outliers5)

  cl <- ClusterFinder1(data5,sigma=0.4,radius=1) 
  plot(data5,col=as.numeric(cl$types),pch=as.numeric(cl$types))
  legend(1,1,legend=levels(cl$types),xjust=1,col=1:length(levels(cl$types)),
             pch=1:length(levels(cl$types)))

}
\keyword{multivariate}
