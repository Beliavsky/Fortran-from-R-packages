\name{endmemberCoordinates}
\alias{endmemberCoordinates}
\alias{endmemberCoordinatesInv}
\alias{endmemberCoordinates.default}
\alias{endmemberCoordinates.acomp}
\alias{endmemberCoordinates.aplus}
\alias{endmemberCoordinates.rplus}
\alias{endmemberCoordinatesInv}
\alias{endmemberCoordinatesInv.rmult}
\alias{endmemberCoordinatesInv.acomp}
\alias{endmemberCoordinatesInv.aplus}
\alias{endmemberCoordinatesInv.rcomp}
\alias{endmemberCoordinatesInv.rplus}
\title{Recast amounts as mixtures of end-members}
\description{
  Computes the convex combination of amounts as mixtures of \code{endmembers}
  to explain \code{X} as well as possible.
}
\usage{
    endmemberCoordinates(X,\dots)
    endmemberCoordinatesInv(K,endmembers,\dots)
    \method{endmemberCoordinates}{default}(X,
                         endmembers=diag(gsi.getD(X)), \dots)
    \method{endmemberCoordinates}{acomp}(X,
                    endmembers=clrInv(diag(gsi.getD(X))),\dots)
    \method{endmemberCoordinates}{aplus}(X,endmembers,\dots)
    \method{endmemberCoordinates}{rplus}(X,endmembers,\dots)
    \method{endmemberCoordinatesInv}{rmult}(K,endmembers,\dots)
    \method{endmemberCoordinatesInv}{acomp}(K,endmembers,\dots)
    \method{endmemberCoordinatesInv}{rcomp}(K,endmembers,\dots)
    \method{endmemberCoordinatesInv}{aplus}(K,endmembers,\dots)
    \method{endmemberCoordinatesInv}{rplus}(K,endmembers,\dots)
          }
\arguments{
  \item{X}{a data set of amounts or compositions, to be represented in
    as convex combination of the \code{endmembers} in the given geometry}
  \item{K}{weights of the \code{endmembers} in the convex combination}
  \item{endmembers}{a dataset of compositions of the same class as X. 
    The number of endmembers given must not exceed the dimension of 
    the space plus one.}
  \item{\dots}{currently unused}
}
\value{
  The \code{endmemberCoordinates} functions give a \code{\link{rmult}}
  data set with the weights (a.k.a. barycentric coordinates) allowing 
  to build \code{X} as good as possible as a convex combination 
  (a mixture) from \code{endmembers}. The result is of class rmult
  because there is no guarantee that the resulting weights are positive
  (although they sum up to one).
  \cr
  The \code{endmemberCoordinatesInv} functions reconstruct the convex
  combination from the weights \code{K} and the given
  \code{endmembers}. The class of \code{endmembers} determines the
  geometry chosen and the class of the result.
}
\details{
  The convex combination is performed in the respective geometry. This
  means that, for rcomp objects, positivity of the result is only guaranteed 
  with endmembers corresponding to extremal individuals of the sample, or 
  completely outside its hull. Note also that, in acomp geometry, the 
  endmembers must necessarily be outside the hull. 
  \cr
  The main idea behind this functions is that the composition actually 
  observed came from a convex combination of some extremal
  compositions, specified by \code{endmembers}. Up to now, this is considered as
  meaningful only in rplus geometry, and under some special circumstances, 
  in rcomp geometry. It is not meaningful in terms of mass conservation 
  in acomp and aplus geometries, because these geometries do not preserve
  mass: whether such an operation has an interpretation is still a matter of 
  debate. In rcomp geometry, the convex combination is dependent on the units of
  measurements, and will be completely different for volume and mass \%. 
  Even more, it is valid only if the whole composition is observed (!). 
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de},
  Raimon Tolosana-Delgado}

\references{
  Shurtz, Robert F., 2003. Compositional geometry and mass conservation.
  Mathematical Geology 35 (8), 972--937.

  
}
\examples{
data(SimulatedAmounts)
ep <- aplus(rbind(c(2,1,2),c(2,2,1),c(1,2,2)))
# mix the endmembers in "ep" with weights given by "sa.lognormals"
dat <- endmemberCoordinatesInv(acomp(sa.lognormals),acomp(ep))
par(mfrow=c(1,2))
plot(dat)
  plot(acomp(ep),add=TRUE,col="red",pch=19)
# compute the barycentric coordinates of the mixture in the "end-member simplex"
plot( acomp(endmemberCoordinates(dat,acomp(ep))))

dat <- endmemberCoordinatesInv(rcomp(sa.lognormals),rcomp(ep))
plot(dat)
plot( rcomp(endmemberCoordinates(dat,rcomp(ep))))

dat <- endmemberCoordinatesInv(aplus(sa.lognormals),aplus(ep))
plot(dat)
plot( endmemberCoordinates(dat,aplus(ep)))

dat <- endmemberCoordinatesInv(rplus(sa.lognormals),rplus(ep))
plot(dat)
plot(endmemberCoordinates(rplus(dat),rplus(ep)))


}
\keyword{multivariate}
