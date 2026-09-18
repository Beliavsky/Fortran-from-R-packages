\name{kerndensmeasure}

\alias{kerndensmeasure}


\title{Statistic measuring closeness to symmetric unimodal distribution}

\description{
  Density-based distance between one-dimensional data and a unimodal
  symmetric distribution about the data mean based on Pons (2013, p.79),
  adapted by Hennig and Coretto (2021), see details.
} 

\usage{
  kerndensmeasure(x,weights=rep(1,nrow(as.matrix(x))),maxq=qnorm(0.9995),
                  kernn=100)
}


\arguments{
  \item{x}{vector. One-dimensional dataset.} 
  \item{weights}{non-negative vector. Relative weights of observations
    (will be standardised to sup up to one internally).}
  \item{maxq}{densities are evaluated between \code{mean(x)-maxq} and
    \code{mean(x)+maxq}.}
  \item{kernn}{integer. Number of points at which the density is
    evaluated, input parameter \code{n} of \code{\link{density}}. This
  should be even.}
}
   
\details{
  Function \code{\link{density}} is used in order to compute a kernel
  density estimator from the data. The \verb{kernn} values of the
  density are then ordered from the pargest to the smallest. Beginning
  from the largest to the smallest, pairs of two values are formed
  (largest and largest biggest, third and fourth largest, and so
  on). Each pair is replaced by two copies of the average of the two
  values. Then on each side of the mean one of each copy is placed from
  the biggest to the smallest, and this produces a symmetric density
  about the mean. The the root mean squared difference between this and
  the original density is computed.
}
  
\value{
  A list with components \code{cp, cpx, measure}.

  \item{cp}{vector of generated symmetric density values from largest to
    smallest (just one copy, sp \code{kernn/2} values).}
  \item{cpx}{\code{y}-component of \code{\link{density}}-output.}
  \item{measure}{root mean squared difference between the densities.}
}

\references{
  Hennig, C. and P.Coretto (2021). An adequacy approach for deciding the
  number of clusters for OTRIMLE robust Gaussian mixture based
  clustering. To appear in \emph{Australian and New Zealand Journal of
  Statistics}, \url{https://arxiv.org/abs/2009.00921}.

  Pons, O. (2013). \emph{Statistical Tests of Nonparametric Hypotheses:
  Asymptotic Theory}. World Scientific, Singapore.
}


\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}}

   
\seealso{
  \code{\link{kerndensp}} 
}

\examples{
  set.seed(124578)
  x <- runif(20)
  str(kerndensmeasure(x))
}

\keyword{htest}


