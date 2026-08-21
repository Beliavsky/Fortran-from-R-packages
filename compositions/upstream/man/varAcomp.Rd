\name{var.acomp}
\alias{var}
\alias{var.default}
\alias{var.acomp}
\alias{var.rcomp}
\alias{var.aplus}
\alias{var.rplus}
\alias{var.rmult}
\alias{cov}
\alias{cov.default}
\alias{cov.acomp}
\alias{cov.rcomp}
\alias{cov.aplus}
\alias{cov.rplus}
\alias{cov.rmult}
\title{Variances and covariances of amounts and compositions}
\description{
  Compute the (co)variance matrix in the several approaches of compositional
 and amount data analysis.
}
\usage{
  var(x,\dots)
    \method{var}{default}(x, y=NULL, na.rm=FALSE, use, \dots)
    \method{var}{acomp}(x,y=NULL,\dots,robust=getOption("robust"), 
       use="all.obs",giveCenter=FALSE)
    \method{var}{rcomp}(x,y=NULL,\dots,robust=getOption("robust"), 
       use="all.obs",giveCenter=FALSE)
    \method{var}{aplus}(x,y=NULL,\dots,robust=getOption("robust"), 
       use="all.obs",giveCenter=FALSE)
    \method{var}{rplus}(x,y=NULL,\dots,robust=getOption("robust"), 
       use="all.obs",giveCenter=FALSE)
    \method{var}{rmult}(x,y=NULL,\dots,robust=getOption("robust"), 
       use="all.obs",giveCenter=FALSE)
  cov(x,y=x,\dots)
    \method{cov}{default}(x, y=NULL, use="everything", 
       method=c("pearson", "kendall", "spearman"), \dots)
    \method{cov}{acomp}(x,y=NULL,\dots,robust=getOption("robust"), 
       use="all.obs",giveCenter=FALSE)
    \method{cov}{rcomp}(x,y=NULL,\dots,robust=getOption("robust"), 
       use="all.obs",giveCenter=FALSE)
    \method{cov}{aplus}(x,y=NULL,\dots,robust=getOption("robust"), 
       use="all.obs",giveCenter=FALSE)
    \method{cov}{rplus}(x,y=NULL,\dots,robust=getOption("robust"), 
       use="all.obs",giveCenter=FALSE)
    \method{cov}{rmult}(x,y=NULL,\dots,robust=getOption("robust"), 
    use="all.obs",giveCenter=FALSE)
          }
\arguments{
  \item{x}{a dataset, eventually of amounts or compositions}
  \item{y}{a second dataset, eventually of amounts or compositions}
  \item{na.rm}{see \code{stats::var}}%{see \code{\link[stats]{var}}}
  \item{use}{see \code{stats::var}}%{see \code{\link[stats]{var}}}
  \item{method}{see \code{stats::cov}}%{see \code{\link[stats]{cov}}}
  \item{\dots}{further arguments to \code{stats::var} %\code{\link[stats]{var}}
    e.g. \code{use}}
  \item{robust}{
    A description of a robust estimator. FALSE for the classical
    estimators. See \link{robustnessInCompositions} for
    further details.
  }
  \item{giveCenter}{If TRUE the center used in the variance calculation
    is reported as a "center" attribute. This is especially necessary
    for robust estimations, where a reasonable center can not be
    computed independently for the me variance calculation.}

}
\value{
  The variance matrix of x or the covariance matrix of x and
  y.}
\details{
  The basic functions of %\code{\link[stats]{var}},\code{\link[stats]{cov}} 
  \code{stats::var} and \code{stats::cov}
  are turned to
  S3-generics. The original versions are copied to the default
  method. This allows us to introduce generic methods to handle
  variances and covariances of other data types, such as amounts or
  compositions.
  \cr
  If classed amounts or compositions are involved, they are transformed
  with their corresponding transforms, using the centered default
  transform (\code{\link{cdt}}). That implies that the variances have to
  be interpreded in a log scale level for  \code{\link{acomp}} and
  \code{\link{aplus}}.
  \cr
  We should be aware that variance matrices of compositions 
  (\code{\link{acomp}} and \code{\link{rcomp}}) are
  singular. They can be transformed to the correponding nonsingular
  variances of ilr or ipt-space by \code{\link{clrvar2ilr}}.

  In R versions older than v2.0.0, 
  %\code{\link[stats]{var}} and \code{\link[stats]{cov}}
  \code{stats::var} and \code{stats::cov}
  were defined in package ``base'' instead of in ``stats''.
  This might produce some misfunction.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{cdt}}, \code{\link{clrvar2ilr}}, \code{\link{clo}},
\code{\link{mean.acomp}}, \code{\link{acomp}}, \code{\link{rcomp}},
\code{\link{aplus}}, \code{\link{rplus}}, \code{\link{variation}}
}

\examples{
data(SimulatedAmounts)
meanCol(sa.lognormals)
var(acomp(sa.lognormals))
var(rcomp(sa.lognormals))
var(aplus(sa.lognormals))
var(rplus(sa.lognormals))
cov(acomp(sa.lognormals5[,1:3]),acomp(sa.lognormals5[,4:5]))
cov(rcomp(sa.lognormals5[,1:3]),rcomp(sa.lognormals5[,4:5]))
cov(aplus(sa.lognormals5[,1:3]),aplus(sa.lognormals5[,4:5]))
cov(rplus(sa.lognormals5[,1:3]),rplus(sa.lognormals5[,4:5]))
cov(acomp(sa.lognormals5[,1:3]),aplus(sa.lognormals5[,4:5]))

svd(var(acomp(sa.lognormals)))

}
\keyword{multivariate}
