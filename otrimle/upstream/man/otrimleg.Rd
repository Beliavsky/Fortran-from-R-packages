\name{otrimleg}
\alias{otrimleg}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{OTRIMLE for a range of numbers of clusters with density-based cluster quality statistic}
\description{
  Computes Optimally Tuned Robust Improper Maximum Likelihood Clustering
  (OTRIMLE), see \code{\link{otrimle}},
  together with the
  density-based cluster quality statistics Q (Hennig and Coretto 2021)
  for a range of values of the number of clusters. 
}

\usage{
otrimleg(dataset, G=1:6, multicore=TRUE, ncores=detectCores(logical=FALSE)-1,
   erc=20, beta0=0, fixlogicd=NULL, monitor=1, dmaxq=qnorm(0.9995))
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dataset}{something that can be coerced into an
    observations times variables matrix. The dataset.}
  \item{G}{vector of integers (normally starting from 1). Numbers of
    clusters to be considered.}
  \item{multicore}{logical. If \code{TRUE}, parallel computing is used
    through the function \code{\link{mclapply}} from package
    \code{parallel}; read warnings there if you intend to use this; it
    won't work on Windows.}
  \item{ncores}{integer. Number of cores for parallelisation.}
  \item{erc}{A number larger or equal than one specifying the maximum
    allowed ratio between within-cluster covariance matrix
    eigenvalues. See \code{\link{otrimle}}.}
  \item{beta0}{A non-negative constant, penalty term for noise, to be
    passed as \code{beta} to  \code{\link{otrimle}}, see documentation
    there.}
  \item{fixlogicd}{numeric of \code{NULL}. Value for the
    logarithm of the improper
    constant density \code{logicd}, see
    \code{\link{rimle}}, which is run instead of \code{\link{otrimle}}
    if this is not \code{NULL}. \code{NULL} means that \code{\link{otrimle}}
    determines it from the data.}
  \item{monitor}{0 or 1. If 1, progress messages are printed on screen.}
  \item{dmaxq}{numeric. Passed as \code{maxq} to
    \code{\link{kerndensmeasure}}. The interval considered for the
    one-dimensional density estimator is \code{(-maxq,maxq)}.}
}

\details{
  For estimating the number of clusters this is meant to be called by
  \code{\link{otrimlesimg}}. The output of \code{\link{otrimleg}} is not
  meant to be used directly for estimating the number of clusters, see
  Hennig and Coretto (2021).}  

\value{
  \code{otrimleg} returns a list
  containing the components \code{solution, iloglik, ibic, criterion,
    logicd, noiseprob, denscrit, ddpm}. All of these are lists or
  vectors of which the component number is the number of clusters.
  \item{solution}{list of output objects of \code{\link{otrimle}} or
    \code{\link{rimle}}.}
  \item{iloglik}{vector of improper likelihood values from
    \code{\link{otrimle}}.}
  \item{ibic}{vector of improper BIC-values (small is good) computed
    from \code{iloglik} and the numbers of parameters. Note that the
    behaviour of the improper likelihood is not compatible with the standard
    use of the BIC, so this is experimental and should not be trusted
    for choosing the number of clusters.}
  \item{criterion}{vector of values of OTRIMLE criterion, see
    \code{\link{otrimle}}.}
  \item{noiseprob}{vector of estimated noise proportions,
    \code{exproportion[1]} from \code{\link{otrimle}}.}
  \item{denscrit}{vector of density-based cluster quality statistics Q
    (Hennig and Coretto 2021) as provided by the
    \code{measure}-component of
    \code{\link{kerndensmeasure}}.}
  \item{ddpm}{list of the vector of cluster-wise density-based cluster
    quality measures as provided by the
    \code{ddpm}-component of \code{\link{kerndensmeasure}}.}
}

\references{
  Coretto, P.  and C. Hennig (2016).  Robust improper maximum
  likelihood: tuning, computation, and a comparison with other
  methods for robust Gaussian clustering.  \emph{Journal of the American
  Statistical Association}, Vol. 111(516), pp. 1648-1659.  doi:
  10.1080/01621459.2015.1100996

  P. Coretto and C. Hennig (2017).  Consistency, breakdown
  robustness, and algorithms for robust improper maximum likelihood
  clustering.  \emph{Journal of Machine Learning Research}, Vol. 18(142),
  pp. 1-39.  \url{https://jmlr.org/papers/v18/16-382.html}

  Hennig, C. and P.Coretto (2021). An adequacy approach for deciding the
  number of clusters for OTRIMLE robust Gaussian mixture based
  clustering. To appear in \emph{Australian and New Zealand Journal of
  Statistics}, \url{https://arxiv.org/abs/2009.00921}.
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}}

\seealso{
  \code{\link{otrimle}}, \code{\link{rimle}}, \code{\link{otrimlesimg}},
  \code{\link{kerndensmeasure}} 
}

\examples{
   data(banknote)
   selectdata <- c(1:30,101:110,117:136,160:161)
   x <- banknote[selectdata,5:7]
   obanknote <- otrimleg(x,G=1:2,multicore=FALSE)
}

\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{robust}

