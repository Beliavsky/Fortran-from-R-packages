\name{otrimle}

\alias{otrimle}

\alias{print.otrimle}

\title{Optimally Tuned Robust Improper Maximum Likelihood Clustering}

\description{
  \code{otrimle} searches for \code{G} approximately Gaussian-shaped
  clusters with/without noise/outliers. The method's tuning controlling  the noise
  level is adaptively chosen based on the data.
}


\usage{
otrimle(data, G, initial = NULL, logicd = NULL, npr.max = 0.5, erc = 20,
beta = 0, iter.max = 500, tol = 1e-06, ncores = NULL, monitor = TRUE)
}


\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Rows correspond
    to observations and columns correspond to variables. Categorical
    variables and \code{NA} values are not allowed.
  }
  \item{G}{
    An integer specifying the number of clusters.
  }
  \item{initial}{An integer vector specifying the initial cluster
    assignment  with \code{0} denoting noise/outliers. If \code{NULL} (default)
    initialization is performed using \code{\link{InitClust}}.
  }
  \item{logicd}{
    A vector defining a grid of \emph{log(icd)} values, where \emph{icd}
    denotes the  improper constant density. If \code{logicd=NULL}
    a default grid is considered. A pure Gaussian Mixture Model fit
    (obtained when \code{log(icd)=-Inf}) is  included in the default
    search path.
  }
  \item{npr.max}{
    A number in \code{[0,1)} specifying the maximum proportion of
    noise/outliers. This defines the \emph{noise proportion
      constraint}. If \code{npr.max=0} a single solution without noise
    component is computed (corresponding to \code{logicd = -Inf}.
  }
  \item{erc}{
    A number \code{>=1} specifying the maximum allowed ratio between
    within-cluster covariance matrix eigenvalues. This defines the
    \emph{eigenratio constraint}. \code{erc=1} enforces spherical clusters
    with equal covariance matrices. A large \code{erc} allows for large
    between-cluster covariance discrepancies. In order to facilitate the
    setting of \code{erc}, it is suggested to scale the columns of
    \code{data} (see \code{\link[base]{scale}}) whenever measurement units
    of the different variables are grossly incompatible.
  }
  \item{beta}{
    A non-negative constant. This is the \emph{beta} penalty coefficient
    introduced in Coretto and Hennig (2016). 
  }
  \item{iter.max}{
    An integer value specifying the maximum number of iterations allowed
    in the underlying ECM-algorithm.
  }
  \item{tol}{
    Stopping criterion for the underlying ECM-algorithm. An ECM iteration
    stops if two successive improper log-likelihood values are within
    \code{tol}.
  }
  \item{ncores}{
    an integer value defining the number of cores used for parallel
    computing. When \code{ncores=NULL} (default), the number \code{r} of
    available cores is detected, and \code{(r-1)} of them are used (See
    \emph{Details}). If \code{ncores=1} no parallel backend is started.
  }
  \item{monitor}{logical. If \code{TRUE} progress  messages are printed on
    screen.
  }
}


\details{
  The \code{otrimle} function computes the OTRIMLE solution based on  the
  ECM-algorithm (expectation conditional maximization algorithm)
  proposed in Coretto and Hennig (2017).

  The otrimle criterion is minimized over the \code{logicd} grid of
  \code{log(icd)} values using parallel computing based on the
  \code{\link[foreach]{foreach}}.
  Note that, depending on the BLAS/LAPACK setting, increasing \code{ncores} may not
  produce the desired reduction in computing time.  The latter
  happens when optimized linear algebra routines are in use (e.g.
  OpenBLAS, Intel Math Kernel Library (MKL), etc.). These optimized shared
  libraries  already implement multithreading. Therefore, in this case
  increasing \code{ncores} may only reduce the computing time marginally.

  Occasionally,  there may be datasets for which the function does not provide a
  solution based on default arguments. This corresponds to
  \code{code=0} and \code{flag=1} or \code{flag=2} in the output (see
  \emph{Value}-section below).  This usually happens when some (or all) of the
  following circumstances occur: (i) \code{erc} is too
  large; (ii) \code{npr.max} is too large; (iii) choice of the initial
  partition. Regarding (i) and (ii) it is not possible to give numeric
  references because whether these numbers are too large/small
  strongly depends on the sample size and the dimensionality of the
  data. References given below explain the relationship between
  these quantities. \cr

  It is suggested to try the following whenever a \code{code=0}
  non-solution occurs. Set the \code{logicd} range wide enough
  (e.g. \code{logicd=seq(-500,-5, length=50)}), choose \code{erc=1},
  and a low choice of \code{npr.max} (e.g. \code{npr.max=0.02}).
  Monitor the solution with the criterion profiling plot
  (\code{\link{plot.otrimle}}). According to the criterion profiling
  plot change \code{logicd}, and increase \code{erc} and \code{npr.max}
  up to the point when a "clear" minimum in the criterion profiling plot
  is obtained. If this strategy does not work it is suggested to
  experiment with a different initial partitions (see \code{initial}
  above).

  TBA: Christian may add something about the \code{beta} here.
  
  The \code{pi} object returned by the \code{rimle} function (see
  \emph{Value}) corresponds to the vector of \code{pi} parameters in
  the underlying pseudo-model (1) defined in Coretto and Hennig (2017).
  With \code{logicd = -Inf} the \code{rimle} function approximates the
  MLE for the \emph{plain} Gaussian mixture model with eigenratio
  covariance regularization, in this case the the first element of the
  \code{pi} vector is set to zero because the noise component is not
  considered.  In general, for iid sampling from finite mixture models
  context, these \emph{pi} parameters define expected clusters'
  proportions. Because of the noise proportion constraint in the RIMLE,
  there are situations where this connection may not happen in
  practice. The latter is likely to happen when both \code{logicd} and
  \code{npr.max} are large. Therefore, estimated expected clusters' proportions
  are reported in the \code{exproportion} object of the \code{rimle}
  output, and these are computed based on the improper posterior
  probabilities given in \code{tau}.  See Coretto and Hennig (2017) for
  more discussion on this.
  
  An earlier approximate version of the algorithm was originally
  proposed in Coretto and Hennig (2016). Software for the original
  version of the algorithm can be found in the supplementary materials
  of Coretto and Hennig (2016).
}



\value{
  An S3 object of class \code{'otrimle'} providing the optimal (according to
  the OTRIMLE criterion) clustering. Output components are as follows:

  \item{code}{
    An integer indicator for the convergence.
    \code{code=0} if no solution is found (see \emph{Details});
    \code{code=1} if at the optimal icd value the corresponding EM-algorithm did not
    converge within \code{em.iter.max};
    \code{code=2} convergence is fully achieved.
  }
  \item{flag}{
    A character string containing one or more flags related to
    the EM iteration at the optimal icd.
    \code{flag=1}  if it was not possible to prevent the numerical
    degeneracy of improper posterior probabilities (\code{tau} value
    below).
    \code{flag=2} if enforcement of the \emph{noise proportion constraint}
    failed for numerical reasons.
    \code{flag=3} if the  \emph{noise proportion constraint} has been
    successfully applied at least once.
    \code{flag=4} if the  \emph{eigenratio constraint} has been
    successfully applied at least once.
  }
  \item{iter}{
    Number of iterations performed in the underlying  EM-algorithm at the
    optimal \code{icd}.
  }
  \item{logicd}{
    Resulting value of the optimal \code{log(icd)}.
  }
  \item{iloglik}{
    Resulting value of the improper likelihood.
  }
  \item{criterion}{
    Resulting value of the OTRIMLE criterion.
  }
  \item{pi}{
    Estimated vector of the \code{pi} parameters of the underlying pseudo-model (see \emph{Details}).
  }
  \item{mean}{
    A matrix of dimension \code{ncol(data) x G} containing the
    mean parameters  of each cluster (column-wise).
  }
  \item{cov}{
    An array of size \code{ncol(data) x ncol(data) x G}
    containing the covariance matrices of each cluster.
  }
  \item{tau}{A matrix of dimension \code{nrow(data) x {1+G}} where
    \code{tau[i, 1+j]} is the estimated (improper) posterior probability that
    the \emph{i}th observation belongs to the \emph{j}th cluster.
    \code{tau[i,1]} is the estimated (improper) posterior probability that
    \emph{i}th observation belongs to the noise component.
  }
  \item{smd}{
    A matrix of dimension  \code{nrow(data) x G } where \code{smd[i,j]}
    is the squared Mahalanobis distance of \code{data[i,]} from
    \code{mean[,j]} according to \code{cov[,,j]}.
  }
  \item{cluster}{
    A vector of integers  denoting cluster assignments for each
    observation. It's \code{0} for observations assigned to  noise/outliers.
  }
  \item{size}{
    A vector of integers with sizes (counts) of each cluster.
  }
  \item{exproportion}{
    A vector of estimated expected clusters' proportions (see \emph{Details}).
  }
  \item{optimization}{
    A data.frame with  the OTRIMLE optimization profiling. For each
    value of \code{log(icd)} explored by the algorithm the data.frame
    stores \code{logicd, criterion, iloglik, code, flag}
    (defined above), and \code{enpr} that denotes the expected noise
    proportion.
  }
}







\section{References}{
  Coretto, P.  and C. Hennig (2016).
  Robust improper maximum likelihood: tuning, computation, and a comparison with
  other methods for robust Gaussian clustering.
  \emph{Journal of the American Statistical Association}, Vol. 111(516),
  pp. 1648-1659. 
  \doi{10.1080/01621459.2015.1100996}

  P. Coretto and C. Hennig (2017).
  Consistency, breakdown robustness, and algorithms for robust improper
  maximum  likelihood clustering.
  \emph{Journal of Machine Learning Research}, Vol. 18(142), pp. 1-39. 
  \href{https://jmlr.org/papers/v18/16-382.html}{https://jmlr.org/papers/v18/16-382.html}
}



\author{Pietro Coretto
    \email{pcoretto@unisa.it}
    \url{https://pietro-coretto.github.io}
}




\seealso{
  \code{\link{plot.otrimle}},
  \code{\link{InitClust}},
  \code{\link{rimle}},
}







\examples{
## Load  Swiss banknotes data
data(banknote)
x <- banknote[,-1]

## Perform otrimle clustering with default arguments
set.seed(1)
a <- otrimle(data=x, G=2, logicd=c(-Inf, -50, -10), ncores=1)

## Plot clustering
plot(a, data=x, what="clustering")

## Plot OTRIMLE criterion profiling
plot(a, what="criterion")

## Plot Improper log-likelihood profiling
plot(a, what="iloglik")

## P-P plot of the clusterwise empirical weighted squared Mahalanobis
## distances against the target distribution pchisq(, df=ncol(data))
plot(a, what="fit")
plot(a, what="fit", cluster=1)



## Perform the same otrimle as before with non-zero penalty
set.seed(1)
b <- otrimle(data=x, G=2, beta = 0.5, logicd=c(-Inf, -50, -10), ncores=1)

## Plot clustering
plot(b, data=x, what="clustering")

## Plot OTRIMLE criterion profiling
plot(b, what="criterion")

## Plot Improper log-likelihood profiling
plot(b, what="iloglik")

## P-P plot of the clusterwise empirical weighted squared Mahalanobis
## distances against the target distribution pchisq(, df=ncol(data))
plot(b, what="fit")
plot(b, what="fit", cluster=1)





   
\dontrun{
## Perform the same example using the finer default grid of logicd
## values using multiple cores
##
a <- otrimle(data = x, G = 2)

## Inspect the otrimle criterion-vs-logicd
plot(a, what = 'criterion')

## The minimum occurs at  a$logicd=-9, and exploring a$optimization it
## cane be seen that the interval [-12.5, -4] brackets the optimal
## solution. We search with a finer grid located around the minimum
##
b <- otrimle(data = x, G = 2, logicd = seq(-12.5, -4, length.out = 25))

## Inspect the otrimle criterion-vs-logicd
plot(b, what = 'criterion')

## Check the difference between the two clusterings
table(A = a$cluster, B = b$cluster)

## Check differences in estimated parameters
##
colSums(abs(a$mean - b$mean))               ## L1 distance for mean vectors
apply({a$cov-b$cov}, 3, norm, type = "F")   ## Frobenius distance for covariances
c(Noise=abs(a$npr-b$npr), abs(a$cpr-b$cpr)) ## Absolute difference for proportions
}
}

