\name{ClusterFinder1}
\alias{ClusterFinder1}
\alias{ClusterFinder1.acomp}
\title{Heuristics to find subpopulations of outliers}
\description{
  The ClusterFinder is a heuristic to find subpopulations of outliers
  essentially by looking for secondary modes in a density
  estimate. 
}
\usage{
ClusterFinder1(X,...)
\method{ClusterFinder1}{acomp}(X,...,sigma=0.3,radius=1,asig=1,minGrp=3,
                                 robust=TRUE)
          }
\arguments{
  \item{X}{the dataset to be clustered}
  \item{\dots}{Further arguments to \code{MahalanobisDist(X,...,robust=robust,pairwise=TRUE)} }
  \item{sigma}{numeric: The Bandwidth of the density estimation kernel
    in a robustly Mahalanobis transformed space. (i.e. in the transform,
where the main group has unit variance)}
  \item{radius}{The minimum size of a cluster in a robustly Mahalanobis
transformed space. (i.e. in the transform, where the main group has unit
variance)}
  \item{asig}{a scaling factor for the geometry of the robustly
    Mahalanobis transformed space when computing the likelihood of an
    observation to belong to group (under a Gaussian assumption). Higher
    values }
  \item{minGrp}{the minimum size of group to be used. Smaller groups are
    treated as single outliers}
  \item{robust}{A robustness description for estimating the variance of
    the main group. FALSE is probably not a usefull value. However later
    other robustness techniques than mcd might be usefull. \code{TRUE}
    just picks the default method of the package.
  }
}
\value{
  A list
   \item{types}{a factor representing the group assignments, when the
     small groups are ignored}
   
   \item{typesTbl}{a table giving the number of members in each of these
     groups}
  \item{groups}{a factor representing the found group assignments}
  \item{isMax}{a logical vector indicating for each observation,whether
    it represent a local maximum in the density estimate.}

  \item{prob}{the infered probability to belong to the different groups
    given as an acomp composition.}
  \item{nmembers}{a tabel giving the number of members of each group}
  \item{density}{the density estimated in each observation location}
  \item{likeli}{The infered likelihood see this observation, for each of
    the groups}
}
\details{
  See \link{outliersInCompositions} for a comprehensive introduction
into the outlier
  treatment in compositions.
  \cr
  The ClusterFinder is labeled with a number to make clear that this is
  just an implementation of some heuristic and not based on some eternal
  truth. Other might give better Clusterfinders. 
  \cr
  Unlike other Clustering Algorithms the basic model of this
  algorithm assumes that there is one dominating subpopulation and an
  unkown number of smaller subpopulations with a similar covariance
  structure but a different mean. The algorithm thus first estimates the
  covariance structure of the main population by a robust location scale
  estimator. Then it uses a simplified (Gaussian) kernel density
estimator to find
  nonrandom secondary modes. The it tries to a assign the different
  observations according to discrimination analysis model to the
  different modes. Groups under a given size are considered as single
  outliers forming a seperate group. In this way the number of clusters
  is kept low even if there are many erratic measurements in the dataset.
  \cr
  The main use of the
  clusters is descriptive plotting. The advantage of these cluster
  against other cluster techniques like k-mean or hclust is that it does
  not tear appart the central mass of the data, as these methods do to
  make the clusters as compact as possible.

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{hclust}}, \code{\link{kmeans}}
}

\examples{
data(SimulatedAmounts)
  cl <- ClusterFinder1(sa.outliers5,sigma=0.4,radius=1) 
  plot(sa.outliers5,col=as.numeric(cl$types),pch=as.numeric(cl$types))
  legend(1,1,legend=levels(cl$types),xjust=1,col=1:length(levels(cl$types)),
                     pch=1:length(levels(cl$types)))

}
\keyword{cluster}
