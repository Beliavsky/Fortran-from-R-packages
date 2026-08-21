\name{mvar}
\alias{mvar}
\alias{mvar.default}
\alias{mcov}
\alias{mcov.default}
\alias{mcor}
\alias{mcor.default}
\alias{msd}
\alias{msd.default}
\title{Metric summary statistics of real, amount or compositional data}
\description{
  Compute the metric variance, covariance, correlation or standard deviation.
}
\usage{
mvar(x,\dots)
mcov(x,\dots)
mcor(x,\dots)
msd(x,\dots)
\method{mvar}{default}(x,y=NULL,\dots)
\method{mcov}{default}(x,y=x,\dots)
\method{mcor}{default}(x,y,\dots)
\method{msd}{default}(x,y=NULL,\dots)
}
\arguments{
  \item{x}{a dataset, eventually of amounts or compositions}
  \item{y}{a second dataset, eventually of amounts or compositions}
  \item{\dots}{further arguments to %\code{\link[stats]{var}} or \code{\link[stats]{cov}}.
  \code{stats::var} or \code{stats::cov}.
     Typically a \code{robust=TRUE} argument. 
     e.g. \code{use}}
}

\value{
  a scalar number, informing of the degree of variation/covariation of one/two datasets.
}
\details{
  The metric variance (\code{mvar}) is defined by the trace of the
  variance in the natural geometry of the data, or also by the generalized
  variance in natural geometry. The natural geometry is equivalently
  given by the \code{\link{cdt}} or \code{\link{idt}} transforms.\cr

  The metric standard deviation (\code{msd}) is not the square root
  of the metric variance, but the square root of the mean of the eigenvalues of the 
  variance matrix. In this way it can be interpreted in units of the original
  natural geometry, as the radius of a sperical ball around
  the mean with the same volume as the 1-sigma ellipsoid of the data set.
  \cr
  
  The metric covariance (\code{mvar}) is the sum over the absolute
  singular values of the covariance of two datasets in their respective 
  geometries. It is always positive. The metric covariance of a dataset 
  with itself is its metric variance. The interpretation of a metric 
  covariance is quite difficult, but useful in regression problems.\cr

  The metric correlation (\code{mcor}) is the metric covariance of the
  datasets in their natural geometry normalized to unit variance matrix. It is a
  number between 0 and the smaller dimension of both natural spaces. A
  number of 1 means perfect correlation in 1 dimension, but only partial
  correlations in higher dimensions.
  \cr
  
}
\references{
  Daunis-i-Estadella, J., J.J. Egozcue, and V. Pawlowsky-Glahn
  (2002) Least squares regression in the Simplex on the simplex, \emph{Terra
    Nostra}, Schriften der Alfred Wegener-Stiftung, 03/2003\cr

  Pawlowsky-Glahn, V. and J.J. Egozcue (2001) Geometric approach to
  statistical analysis on the simplex. \emph{SERRA} \bold{15}(5), 384-398\cr
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}
\seealso{
\code{\link{var}}, \code{\link{cov}},
\code{\link{mean.acomp}}, \code{\link{acomp}}, \code{\link{rcomp}},
\code{\link{aplus}}, \code{\link{rplus}} }
\examples{
data(SimulatedAmounts)
mvar(acomp(sa.lognormals))
mvar(rcomp(sa.lognormals))
mvar(aplus(sa.lognormals))
mvar(rplus(sa.lognormals))

msd(acomp(sa.lognormals))
msd(rcomp(sa.lognormals))
msd(aplus(sa.lognormals))
msd(rplus(sa.lognormals))

mcov(acomp(sa.lognormals5[,1:3]),acomp(sa.lognormals5[,4:5]))
mcor(acomp(sa.lognormals5[,1:3]),acomp(sa.lognormals5[,4:5]))
mcov(rcomp(sa.lognormals5[,1:3]),rcomp(sa.lognormals5[,4:5]))
mcor(rcomp(sa.lognormals5[,1:3]),rcomp(sa.lognormals5[,4:5]))

mcov(aplus(sa.lognormals5[,1:3]),aplus(sa.lognormals5[,4:5]))
mcor(aplus(sa.lognormals5[,1:3]),aplus(sa.lognormals5[,4:5]))
mcov(rplus(sa.lognormals5[,1:3]),rplus(sa.lognormals5[,4:5]))
mcor(rplus(sa.lognormals5[,1:3]),rplus(sa.lognormals5[,4:5]))

mcov(acomp(sa.lognormals5[,1:3]),aplus(sa.lognormals5[,4:5]))
mcor(acomp(sa.lognormals5[,1:3]),aplus(sa.lognormals5[,4:5]))
}
\keyword{multivariate}
