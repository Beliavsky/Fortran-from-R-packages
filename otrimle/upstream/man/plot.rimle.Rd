\name{plot.rimle}
\alias{plot.rimle}

\title{
  Plot Methods for RIMLE Objects
}

\description{
  Plot robust model-based clustering results: scatter plot with
  clustering information and cluster fit.
}

\usage{
\method{plot}{rimle}(x, what=c("fit", "clustering"),
     data=NULL, margins=NULL, cluster=NULL, \dots)
}

\arguments{
  \item{x}{
    Output from \code{\link{rimle}}
  }
  \item{what}{
    The type of graph. It can be one of the following:
     \code{"fit"} (default), \code{"clustering"}. See \emph{Details}.
  }
  \item{data}{The data vector, matrix or data.frame (or some
    transformation of them),  used for obtaining the
    \code{'rimle'} object. This is only relevant if
    \code{what="clustering"}.
  }
  \item{margins}{
    A vector of integers denoting the variables (numbers of columns of
  \code{data}) to be used for a \code{pairs}-plot if
    \code{what="clustering"}. When \code{margins=NULL} it is
    set to  \code{1:ncol(data)} (default).
  }
  \item{cluster}{An integer denoting the cluster for which the \emph{fit}
    plot is returned. This is only relevant if \code{what="fit"}.
  }
  \item{\dots}{
    further arguments passed to or from other methods.
  }
}


\value{
  \describe{
   \item{If \code{what="fit"}}{
      The P\code{-}P plot (probability\code{-}probability plot) of the weighted empirical
      distribution function of the Mahalanobis distances of observations
      from clusters' centers against the target distribution. The target
      distribution is the Chi-square distribution with degrees of
      freedom equal to \code{ncol(data)}. The weights are given by the improper posterior
      probabilities. If \code{cluster=NULL} P\code{-}P plots are produced for
      all clusters, otherwise \code{cluster} selects a single P\code{-}P
      plot at times.
    }
    \item{If \code{what="clustering"}}{
      A pairwise scatterplot with cluster memberships. Points
      assigned to the noise/outliers component are denoted by
      \code{'+'}.
    }
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
  \code{\link{otrimle}}
}



\examples{
## Load  Swiss banknotes data
data(banknote)
x <- banknote[,-1]

## Perform rimle clustering with default arguments
set.seed(1)
a <- rimle(data = x, G = 2)
print(a)

## Plot clustering
plot(a, data = x, what = "clustering")

## Plot clustering on selected margins
plot(a, data = x, what = "clustering", margins = 4:6)

## Plot clustering on the first two principal components
z <- scale(x) \%*\%   eigen(cor(x), symmetric = TRUE)$vectors
colnames(z) <- paste("PC", 1:ncol(z), sep = "")
plot(a, data = z, what = "clustering", margins = 1:2)

## Fit plot for all clusters
plot(a, what = "fit")

## Fit plot for cluster 1
plot(a, what = "fit", cluster = 1)
}
