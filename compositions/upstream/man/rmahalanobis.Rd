\name{rMahalanobis}
\alias{pMaxMahalanobis}
\alias{qMaxMahalanobis}
\alias{rMaxMahalanobis}
\alias{pEmpiricalMahalanobis}
\alias{qEmpiricalMahalanobis}
\alias{rEmpiricalMahalanobis}
\alias{pPortionMahalanobis}
\alias{qPortionMahalanobis}
\alias{rPortionMahalanobis}
\alias{pQuantileMahalanobis}
\alias{gsi.pStore}
\title{Compute distributions of empirical Mahalanobis distances based on
simulations}
\description{
  Decissions about outliers are often made based on Mahalanobis
  distances with respect to robustly estimated variances. These function
  deliver the necessary distributions.
}
\usage{
rEmpiricalMahalanobis(n,N,d,...,sorted=FALSE,pow=1,robust=TRUE)
pEmpiricalMahalanobis(q,N,d,...,pow=1,replicates=100,resample=FALSE,
                        robust=TRUE)
qEmpiricalMahalanobis(p,N,d,...,pow=1,replicates=100,resample=FALSE,
                        robust=TRUE)
rMaxMahalanobis(n,N,d,...,pow=1,robust=TRUE)
pMaxMahalanobis(q,N,d,...,pow=1,replicates=998,resample=FALSE,
                        robust=TRUE)
qMaxMahalanobis(p,N,d,...,pow=1,replicates=998,resample=FALSE,
                        robust=TRUE)
rPortionMahalanobis(n,N,d,cut,...,pow=1,robust=TRUE)
pPortionMahalanobis(q,N,d,cut,...,replicates=1000,resample=FALSE,pow=1,
                        robust=TRUE)
qPortionMahalanobis(p,N,d,cut,...,replicates=1000,resample=FALSE,pow=1,
                        robust=TRUE)
pQuantileMahalanobis(q,N,d,p,...,replicates=1000,resample=FALSE,
                        ulimit=TRUE,pow=1,robust=TRUE)
}
\arguments{
  \item{n}{Number of simulations to do.}
  \item{q}{A vector giving quantiles of the distribution}
  \item{p}{A vector giving probabilities. (only a single probility for 
             \code{pQuantileMahalanobis})}
  \item{N}{Number of cases in the dataset.}
  \item{d}{degrees of freedom (i.e. dimension) of the dataset.}
  \item{cut}{A cutting limit. The random variable is the portion of
    Mahalanobis distances lower equal to the cutting limit.}
  \item{\ldots}{further arguments passed to \code{\link{MahalanobisDist}}}
  \item{pow}{the power of the Mahalanobis distance to be used. Higher
    powers can be used to stretch the outlierregion visually.}
  \item{robust}{logical or a robust method description (see
    \code{\link{robustnessInCompositions}}) specifiying how the center
and covariance
    matrix are estimated,if not given.}
  \item{sorted}{
    Specifies a transformation to be applied to the whole sequence of
    Mahalanobis distances: FALSE is no transformation, TRUE sorts the
    entries in ascending order, a numeric vector picks the given entries
    from the entries sorted in ascending order; alternatively a function
    such as \code{max} can be given to directly transform the data.
  }
  \item{replicates}{the number of datasets in the
    Monte-Carlo-Computations used in these routines.}
  \item{resample}{a logical forcing a resampling of the
    Monte-Carlo-Sampling. See details. }
  \item{ulimit}{logical: is this an upper limit of a joint confidence
    bound or a lower limit.}
}
\value{
  The r* functions deliver a vector (or a matrix of row-vectors) of
  simulated value of the given distributions. A total of n values (or
  row vectors) is returned.
  \cr
  The p* functions deliver a vector (of the same length as x) of
  probabilities for random variable of the given distribution to be
  under the given quantil values q.
  \cr
  The q* functions deliver a vector of quantiles corresponding to the
  length of the vector p providing the probabilities.
  
}
\details{
  All the distribution correspond to the distribution under the
  Null-Hypothesis of multivariate joint Gaussian distribution of the
dataset.
  
  The set of empirically estimated Mahalanobis distances of a dataset is
  in the first step a random vector with exchangable but dependent
  entries. The distribution of this vector is given by the
  \code{rEmpiricalMahalanobis} if no sorted argument is given. Please be
  advised that this is not a fixed distribution in a mathematical sense,
  but an implementation dependent distribution incorporating the
  performance of underlying robust spread estimator. As long as no
  sorted argument is given \code{pEmpiricalMahalanobis} and
  \code{qEmpiricalMahalanobis} represent the distribution function and
  the quantile function of a randomly picked element of this
  vector.
  \cr
  If a sorted attribute is given, it specifies a transformation is
  applied to each of the vector prior to processing. Three important
  special cases
  are provided by seperate functions. The MaxMahalanobis functions
  correspond to picking only the larges value. The PortionMahalanobis
  functions correspond to reporting the portion of Mahalanobis distances
  over a cutoff. The QuantileMahalanobis distribution correponds to the
  distribution of the p-quantile of the dataset. 
  \cr
  The Monte-Carlo-Simulations of these
    distributions are rather slow, since for each datum we need to
    simulate a whole dataset and to apply a robust covariance estimator
    to it, which typically itself involves
    Monte-Carlo-Algorithms. Therefore each type of simulations is only
    done the first time needed and stored for later use in the
    environment \code{\link{gsi.pStore}}. With the resampling argument a
    resampling of the cashed dataset can be forced.
  
}
\note{
  Unlike the \code{\link{mahalanobis}} function this function does not
  be default compute the square of the mahalanobis distance. The pow
  option is provided if the square is needed.
  \cr
  The package \pkg{robustbase} is required for using the
  robust estimations.}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{dist}}, \code{\link{OutlierClassifier1}}
}

\examples{
rEmpiricalMahalanobis(10,25,2,sorted=TRUE,pow=1,robust=TRUE)
pEmpiricalMahalanobis(qchisq(0.95,df=10),11,1,pow=2,replicates=1000)
(xx<-pMaxMahalanobis(qchisq(0.95,df=10),11,1,pow=2))
qEmpiricalMahalanobis(0.95,11,2)
rMaxMahalanobis(10,25,4)
qMaxMahalanobis(xx,11,1)
}
\keyword{distribution}
