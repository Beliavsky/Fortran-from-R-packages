
\name{rimle}

\alias{rimle}

\alias{print.rimle}

\title{Robust Improper Maximum Likelihood Clustering}

\description{
  \code{rimle} searches for \code{G} approximately Gaussian-shaped
  clusters with/without noise/outliers. The method's tuning controlling
  the noise level is fixed and is to be provided by the user or will be guessed by
  the function in a rather quick and dirty way (\code{\link{otrimle}}
  performs  a more sophisticated data-driven choice).
}


\usage{
rimle(data, G, initial=NULL, logicd=NULL, npr.max=0.5, erc=20, iter.max=500, tol=1e-6)
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
    A number \code{log(icd)}, where \code{0 <= icd < Inf}
    is the value of the improper constant density (\code{icd}). This is
    the RIMLE's tuning for controlling the size of the noise.
    If \code{logicd=NULL} (default), an \code{icd} value is guessed
    based on the data. A pure Gaussian Mixture Model fit is obtained
    with \code{logicd = -Inf}.
  }
  \item{npr.max}{
    A number in \code{[0,1)} specifying the maximum proportion of
    noise/outliers. This defines the \emph{noise proportion
    constraint}. If \code{npr.max=0} a solution without noise
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
  \item{iter.max}{
    An integer value specifying the maximum number of iterations allowed
    in the  ECM-algorithm (see Details).
  }
  \item{tol}{
    Stopping criterion for the underlying ECM-algorithm. An ECM iteration
    stops if two successive improper log-likelihood values are within
    \code{tol}.
  }
}


\details{
  The \code{rimle} function computes the RIMLE solution using the
  ECM-algorithm proposed in Coretto and Hennig (2017).
  
  There may be datasets for which the function does not provide a
  solution based on default arguments. This corresponds to
  \code{code=0} and \code{flag=1} or \code{flag=2}  in the output (see
  \emph{Value}-section below).  This usually happens when some (or all) of the
  following circumstances occur: (i) \code{log(icd)} is too
  large; (ii) \code{erc} is too large; (iii) npr.max is too large;
  (iv) choice of the initial partition. In these cases it is suggested
  to find a suitable interval of \code{icd} values by using the
  \code{\link{otrimle}} function. The \emph{Details} section of
  \code{\link{otrimle}} suggests several actions to take
  whenever a \code{code=0} non-solution occurs.

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
  \code{npr.max} are large. Therefore, estimated expected clusters'
  proportions are reported in the \code{exproportion} object of the
  \code{rimle} output, and these are computed based on the
  improper posterior probabilities given in \code{tau}.
  See Coretto and Hennig (2017) for more discussion on this. 

  An earlier approximate version of the algorithm was originally
  proposed in Coretto and Hennig (2016). Software for the original
  version of the algorithm can be found in the supplementary materials
  of  Coretto and Hennig (2016).
}



\value{
  An S3 object of class \code{'rimle'}. Output components are as follows:

  \item{code}{
    An integer indicator for the convergence.
    \code{code=0} if no solution is found (see \emph{Details});
    \code{code=1} if the  EM-algorithm did not converge within \code{em.iter.max};
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
    \code{flag=3} if enforcement of the \emph{eigenratio constraint} failed for
    numerical reasons.
    \code{flag=4} if the  \emph{noise proportion constraint} has been
    successfully applied at least once.
    \code{flag=5} if the  \emph{eigenratio constraint} has been
    successfully applied at least once.
  }
  \item{iter}{
    Number of iterations performed in the underlying  EM-algorithm.
  }
  \item{logicd}{
    Value of the \code{log(icd)}.
  }
  \item{iloglik}{
    Value of the improper likelihood.
  }
  \item{criterion}{
    Value of the OTRIMLE criterion.
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
}







\references{
  Coretto, P.  and C. Hennig (2016).
  Robust improper maximum likelihood: tuning, computation, and a comparison with
  other methods for robust Gaussian clustering.
  \emph{Journal of the American Statistical Association}, Vol. 111(516),
  pp. 1648-1659. \doi{10.1080/01621459.2015.1100996}


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
  \code{\link{plot.rimle}},
  \code{\link{InitClust}},
  \code{\link{otrimle}},
}







\examples{
## Load  Swiss banknotes data
data(banknote)
x <- banknote[,-1]

## -----------------------------------------------------------------------------
## EXAMPLE 1:
## Perform RIMLE with default inputs
## -----------------------------------------------------------------------------
set.seed(1)
a <- rimle(data = x, G = 2)
print(a)

## Plot clustering
plot(a, data = x, what = "clustering")

## P-P plot of the clusterwise empirical weighted squared Mahalanobis
## distances against the target distribution pchisq(, df=ncol(data))
plot(a, what = "fit")
plot(a, what = "fit", cluster = 1)



## -----------------------------------------------------------------------------
## EXAMPLE 2:
## Compare solutions for different choices of logicd
## -----------------------------------------------------------------------------
set.seed(1)

## Case 1: noiseless solution, that is fit a pure Gaussian Mixture Model
b1 <- rimle(data = x, G = 2, logicd = -Inf)
plot(b1, data=x, what="clustering")
plot(b1, what="fit")

## Case 2: low noise level
b2 <- rimle(data = x, G = 2, logicd = -100)
plot(b2, data=x, what="clustering")
plot(b2, what="fit")

## Case 3: medium noise level
b3 <- rimle(data = x, G = 2, logicd = -10)
plot(b3, data=x, what="clustering")
plot(b3, what="fit")

## Case 3: large noise level
b3 <- rimle(data = x, G = 2, logicd = 5)
plot(b3, data=x, what="clustering")
plot(b3, what="fit")

}

