\name{logratioVariogram}
\alias{logratioVariogram}
\alias{logratioVariogram.acomp}
\title{Empirical variograms for compositions}
\description{
  Computes the matrix of logratio variograms.
}
\usage{
  logratioVariogram(data, \dots)
  
  \method{logratioVariogram}{acomp}(data,
                          loc,
                          maxdist=max(dist(loc))/2,
                          nbins=20,
                          dists=seq(0,maxdist,length.out=nbins+1),
                          bins=cbind(dists[-length(dists)],dists[-1]),
                          azimuth=0,
                          azimuth.tol=180,
                          comp=data,
                          \dots
                          )
          }
\arguments{
  \item{data}{an acomp compositional dataset}
  \item{\dots}{arguments for generic functionality}
  \item{loc}{a matrix or dataframe providing the observation locations
    of the compositions. Any number of dimension >= 2 is supported.}
  \item{maxdist}{the maximum distance to compute the variogram for.}
  \item{nbins}{The number of distance bins to compute the variogram for}
  \item{dists}{The distances seperating the bins}
  \item{bins}{a matrix with lower and upper limit for the distances of
    each bin. A pair is counted if min<h<=max. min and max are provided
    as columns. bins is computed from maxdist,nbins and dists. If it is
    provided, it is used directly.}
  \item{azimuth}{ For directional variograms the direction, either as an
    azimuth  angle (i.e. a single real number) for 2D
    datasets or a unit vector pointing of the same dimension as the
    locations. The angle is clockwise from North in degree. }
  \item{azimuth.tol}{ The angular tolerance it should be below 90 if a
    directional variogram is intended. }
  \item{comp}{do not use, only provided for backwards compatibility. Use \code{data} instead}
}
\value{
  A list of class \code{"logratioVariogram"}.
  \item{vg}{A nbins x D x D array containing the logratio variograms}
  \item{h}{A nbins x D x D array containing the mean distance the
    value is computed on. }
  \item{n}{A nbins x D x D array containing the
    number of nonmissing pairs used for the corresponding value.}
}
\details{
  The logratio-variogram is the set of variograms of each of the pairwise
  logratios. It can be proven that it carries the same information as a
  usual multivariate variogram. The great advantage is that all the
  funcitions have a direct interpreation and can be estimated even with
  (MAR) missings in the dataset. 
}
\references{
  Tolosana, van den Boogaart, Pawlowsky-Glahn (2009) Estimating and
  modeling variograms of compositional data with occasional missing
  variables in R, StatGis09
  
  Pawlowsky-Glahn, Vera and Olea, Ricardo A. (2004) Geostatistical
  Analysis of Compositional Data, Oxford University Press, Studies in
  Mathematical Geology
  
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{vgram2lrvgram}},
  \code{\link{CompLinModCoReg}},
  \code{\link{vgmFit}}
}
\examples{
\dontrun{
data(juraset)
X <- with(juraset,cbind(X,Y))
comp <- acomp(juraset,c("Cd","Cu","Pb","Co","Cr"))
lrv <- logratioVariogram(comp,X,maxdist=1,nbins=10)
plot(lrv)
}
}
\keyword{multivariate}
