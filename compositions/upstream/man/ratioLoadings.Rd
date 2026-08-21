\name{ratioLoadings}
\alias{relativeLoadings}
\alias{relativeLoadings.princomp.acomp}
\alias{relativeLoadings.princomp.aplus}
\alias{relativeLoadings.princomp.rcomp}
\alias{relativeLoadings.princomp.rplus}
\alias{print.relativeLoadings.princomp.acomp}
\alias{print.relativeLoadings.princomp.aplus}
\alias{print.relativeLoadings.princomp.rcomp}
\alias{print.relativeLoadings.princomp.rplus}
\alias{plot.relativeLoadings.princomp.acomp}
\alias{plot.relativeLoadings.princomp.aplus}
\alias{plot.relativeLoadings.princomp.rcomp}
\alias{plot.relativeLoadings.princomp.rplus}
\title{Loadings of relations of two amounts}
\description{
  In a compositional dataset the relation of two objects can be
  interpreted safer than a single amount. These functions compute,
  display and plot the corresponding pair-information for the various
  principal component analysis results. 
}
\usage{
relativeLoadings(x,...)
\method{relativeLoadings}{princomp.acomp}(x,\dots,log=FALSE,scale.sdev=NA,
                                                  cutoff=0.1)
\method{relativeLoadings}{princomp.aplus}(x,\dots,log=FALSE,scale.sdev=NA,
                                                  cutoff=0.1)
\method{relativeLoadings}{princomp.rcomp}(x,\dots,scale.sdev=NA,
                                                  cutoff=0.1)
\method{relativeLoadings}{princomp.rplus}(x,\dots,scale.sdev=NA,
                                                  cutoff=0.1)
\method{print}{relativeLoadings.princomp.acomp}(x,\dots,cutoff=attr(x,"cutoff"),
                                                  digits=2)
\method{print}{relativeLoadings.princomp.aplus}(x,\dots,cutoff=attr(x,"cutoff"),
                                                 digits=2)
\method{print}{relativeLoadings.princomp.rcomp}(x,\dots,cutoff=attr(x,"cutoff"),
                                                 digits=2)
\method{print}{relativeLoadings.princomp.rplus}(x,\dots,cutoff=attr(x,"cutoff"),
                                                 digits=2)
\method{plot}{relativeLoadings.princomp.acomp}(x,\dots)
\method{plot}{relativeLoadings.princomp.aplus}(x,\dots)
\method{plot}{relativeLoadings.princomp.rcomp}(x,\dots)
\method{plot}{relativeLoadings.princomp.rplus}(x,\dots)
}
\arguments{
  \item{x}{a result from an amount PCA \code{\link{princomp.acomp}}/\code{\link{princomp.aplus}}/\code{\link{princomp.rcomp}}/\code{\link{princomp.rplus}}}
  \item{log}{a logical indicating to use log-ratios instead of ratios}
  \item{scale.sdev}{if not \code{NA}, a number specifying the
    multiple of a standard deviation, used to scale the components}
  \item{cutoff}{a single number. Changes under that (log)-cutoff are not
    displayed}
  \item{digits}{the number of digits to be displayed}
  \item{\dots}{further parameters to internally-called functions}
  }
\value{
  The value is a matrix of type
  \code{"relativeLoadings.princomp.*"}, containing the ratios in the
  compositions represented by the loadings (optionally scaled by the
  standard deviation of the components and \code{scale.sdev}). 
}
\details{
  The relative loadings of components allow a direct interpretation of the effects
  of principal components. For acomp/aplus classes the relation is
  induced by a ratio, which can optionally be log-transformed. For the
  rcomp/rplus-classes the relation is induced by a difference, which
  is meaningless when the units are different.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{princomp.acomp}},   \code{\link{princomp.aplus}},
  \code{\link{princomp.rcomp}},   \code{\link{princomp.rplus}},
  \code{\link[graphics]{barplot}}
}

\examples{
data(SimulatedAmounts)
pc <- princomp(acomp(sa.lognormals5))
pc
summary(pc)
relativeLoadings(pc,log=TRUE)
relativeLoadings(pc)
relativeLoadings(pc,scale.sdev=1)
relativeLoadings(pc,scale.sdev=2)

plot(relativeLoadings(pc,log=TRUE))
plot(relativeLoadings(pc))
plot(relativeLoadings(pc,scale.sdev=1))
plot(relativeLoadings(pc,scale.sdev=2))


}
\keyword{multivariate}
