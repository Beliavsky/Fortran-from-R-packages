\name{ellipses}
\alias{ellipses}
\alias{ellipses.rmult}
\alias{ellipses.acomp}
\alias{ellipses.rcomp}
\alias{ellipses.aplus}
\alias{ellipses.rplus}
\title{Draw ellipses}
\description{
  Draws ellipses from a mean and a variance into a plot.
}
\usage{
  \method{ellipses}{acomp}(mean,var,r=1,\dots,steps=72,
                                         thinRatio=NULL,aspanel=FALSE)
  \method{ellipses}{rcomp}(mean,var,r=1,\dots,steps=72,
                                         thinRatio=NULL,aspanel=FALSE)
  \method{ellipses}{aplus}(mean,var,r=1,\dots,steps=72,thinRatio=NULL)
  \method{ellipses}{rplus}(mean,var,r=1,\dots,steps=72,thinRatio=NULL)
  \method{ellipses}{rmult}(mean,var,r=1,\dots,steps=72,thinRatio=NULL)
          }
\arguments{
  \item{mean}{a compositional dataset or value of means or midpoints of
    the ellipses}
  \item{var}{a variance matrix or a set of variance matrices given by
    \code{var[i,,]} (multiple covariance matrices are not consitently
    implemented as of today). The principal axis of the variance give
    the axis of
  the ellipses, whereas the square-root of the eigenvalues times r give the
  half-diameters of the ellipse.}
  \item{r}{a scaling of the half-diameters}
  \item{\dots}{further graphical parameters}
  \item{steps}{the number of discretisation points to draw the
    ellipses.}
  \item{thinRatio}{The ellipse function now be default plots the whole
    ellipsiod by giving its principle circumferences. However this is
    not reasonable for the thinner directions. If a direction other than
    the first two eigendirections has an eigenvalue not bigger than
    thinRatio*rmax it is not plotted. Thus thinRatio=1 reinstantiates
    the old behavior of the function. Later thinratio=NULL will become the
    default, in which case the projection of the ellipse is
    plotted. However this is not implemented yet. 
  }
  \item{aspanel}{Is the function called as slave to draw in a panel of a
  gsi.pairs plot, or as a user function setting up the plots. }
}
\details{
  The ellipsoid/ellipse drawn is given by the solutions of
  \deqn{(x-mean)^tvar^{-1}(x-mean)=r^2}{(x-mean)^tvar^{-1}(x-mean)=r^2}
  in the respective geometry of the parameter space. Note that
  these ellipses can be added to panel plots (by means of 
  orthogonal projections in the corresponding geometry).

  There are actually three possibilities of drawing a a hyperdimensional
  ellipsoid or ellipse and non of them is perfect.
  \describe{
  \item{thinRatio=1.1}{
    This works like, what was implemented in the older versions of
    compositons, but never correctly documented. It draws an ellipse
    with main axes given by the two largest Eigendirections of the
    \code{var}-Matrix given. 
  }
  \item{thinRatio=0}{
    Draws all the ellipses given by every pair of eigendirections. In
    this way we get a visual impression of the high dimensional
    ellipsoid represend by the variance matrix. However the plots gets
    fastly cluttered in dimensions, when D>4. A 0<thinRatio<1 can avoid
    using eigendirection with small extend (i.e. smaller than
    thinRatio*largest Eigenvalue. 
  }
  \item{thinRatio=NULL}{
    Draws in each Panel a two dimensional ellipse representing the
    marginal variance in the projection of the plot, if var was to be
    interpreted as a variance matrix. This can be seen as some sort of
    projection of the high dimensional ellipsoid, but is not necessarily
    its visual outline. 
  }
  }
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{plot.acomp}}, 
}

\examples{
data(SimulatedAmounts)

plot(acomp(sa.lognormals))
tt<-acomp(sa.lognormals); ellipses(mean(tt),var(tt),r=2,col="red")
tt<-rcomp(sa.lognormals); ellipses(mean(tt),var(tt),r=2,col="blue")

plot(aplus(sa.lognormals[,1:2]))
tt<-aplus(sa.lognormals[,1:2]); ellipses(mean(tt),var(tt),r=2,col="red")
tt<-rplus(sa.lognormals[,1:2]); ellipses(mean(tt),var(tt),r=2,col="blue")

plot(rplus(sa.lognormals[,1:2]))
tt<-aplus(sa.lognormals[,1:2]); ellipses(mean(tt),var(tt),r=2,col="red")
tt<-rplus(sa.lognormals[,1:2]); ellipses(mean(tt),var(tt),r=2,col="blue")
tt<-rmult(sa.lognormals[,1:2]); ellipses(mean(tt),var(tt),r=2,col="green")

}
\keyword{multivariate}
