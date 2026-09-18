\name{kerndensp}
\alias{kerndensp}


\title{Closeness of multivariate distribution to elliptical unimodal distribution}

\description{
  Density- and and principal components-based distance between
  multivariate data and a unimodal
  elliptical distribution about the data mean, see Hennig and Coretto
  (2021). For use in \code{\link{kerndenscluster}}.
} 

\usage{
  kerndensp(x,weights=rep(1,nrow(as.matrix(x))), siglist,maxq=qnorm(0.9995),
            kernn=100)
}



\arguments{
  \item{x}{something that can be coerced into a matrix. Dataset.} 
  \item{weights}{non-negative vector. Relative weights of observations
    (will be standardised to sup up to one internally).}
  \item{siglist}{list with components \code{cov} (covariance matrix),
    \code{center} (mean) and \code{n.obs} (number of observations).}
  \item{maxq}{positive numeric. One-dimensional densities are evaluated
    between \code{mean(x)-maxq} and \code{mean(x)+maxq}.}
  \item{kernn}{integer. Number of points at which the one-dimensional
    density is evaluated, input parameter \code{n} of
    \code{\link{density}}. This should be even.}
}

\details{
  See Hennig and Coretto (2021), Sec. 4.2. \code{\link{kerndensmeasure}}
  is run on the principal components of \code{x}. The resulting measures
  are standardised by \code{\link{kmeanfun}} and \code{\link{ksdfun}}
  and then aggregated as mean square of the positive values, see
  Hennig and Coretto (2021). The PCS is computed by
  \code{\link{princomp}} and will always use \code{siglist} rather than
  statistics computed from \code{x}. 
}
  
\value{
  A list with components \code{cml, cm, pca, stanmeasure, measure}.

  \item{cml}{List of outputs of \code{\link{kerndensmeasure}} for the
    principal components.}
  \item{cm}{vector of \code{measure} components of
    \code{\link{kerndensmeasure}} for the principal components.}
  \item{stanmeasure}{vector of standardised \code{measure} components of
    \code{\link{kerndensmeasure}} for the principal components.}
  \item{pca}{output of \code{\link{princomp}}.} 
  \item{measure}{Final aggregation result.}
}

\references{
  Hennig, C. and P.Coretto (2021). An adequacy approach for deciding the
  number of clusters for OTRIMLE robust Gaussian mixture based
  clustering. To appear in \emph{Australian and New Zealand Journal of
  Statistics}, \url{https://arxiv.org/abs/2009.00921}.
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}

\seealso{
  \code{\link{kerndensmeasure}},  \code{\link{kerndenscluster}}
}

\examples{
  set.seed(124578)
  x <- cbind(runif(20),runif(20))
  siglist <- list(cov=cov(x),center=colMeans(x),n.obs=20) 
  kerndensp(x,siglist=siglist)$measure
}

\keyword{htest}


