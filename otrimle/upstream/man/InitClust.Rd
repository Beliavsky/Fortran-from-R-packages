\name{InitClust}

\alias{InitClust}

\title{Robust Initialization for Model-based Clustering Methods}

\description{
  Computes the initial cluster assignment based on a combination of
  nearest neighbor based noise detection, and agglomerative hierarchical
  clustering based on maximum likelihood criteria for Gaussian
  mixture models.
}

\usage{
 InitClust(data , G , k = 3 , knnd.trim = 0.5 , modelName='VVV')
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
   \item{k}{
     An integer specifying the number of considered nearest neighbors per point
     used for the denoising step (see \emph{Details}).
   }
   \item{knnd.trim}{
     A number in [0,1) which defines the proportion of points
     initialized as noise. Tipically \code{knnd.trim<=0.5} (see \emph{Details}).
   }
   \item{modelName}{
     A character string indicating the covariance model to be used. Possible models are: \cr
     \code{"E"}:   equal variance  (one-dimensional) \cr
     \code{"V"} :  spherical, variable variance (one-dimensional) \cr
     \code{"EII"}: spherical, equal volume \cr
     \code{"VII"}: spherical, unequal volume \cr
     \code{"EEE"}: ellipsoidal, equal volume, shape, and orientation \cr
     \code{"VVV"}: ellipsoidal, varying volume, shape, and orientation (default).\cr
     See \emph{Details}.
   }
 }




 \details{
   The initialization is based on Coretto and Hennig (2017). First, wwo
   steps are performed:\cr

   \emph{Step 1 (denoising step):} for each data point compute its
   \code{k}th\code{-}nearest neighbors
   distance (\code{k-}NND). All points with \code{k-}NND  larger
   than the (1-\code{knnd.trim})\code{-}quantile  of the \code{k-}NND
   are initialized as noise. Intepretaion of
   \code{k} is that:  \code{(k-1)}, but not \code{k}, points close
   together may still be interpreted  as noise or outliers

   \emph{Step 2 (clustering step):} perform the model-based hierarchical
   clustering (MBHC) proposed in Fraley (1998). This step is performed using
   \code{\link[mclust]{hc}}. The input argument \code{modelName} is passed
   to \code{\link[mclust]{hc}}. See \emph{Details} of
   \code{\link[mclust]{hc}} for more details.

   If the previous \emph{Step 2} fails to provide \code{G} clusters each
   containing at least 2 distinct data points, it is replaced with
   classical hirararchical clustering implemented in
   \code{\link[stats]{hclust}}. Finally, if
   \code{\link[stats]{hclust}} fails to provide a valid partition, up
   to ten random partitions are tried.
 }


 \value{An integer vector specifying the initial cluster
    assignment  with \code{0} denoting noise/outliers.
 }


 \section{References}{
   Fraley, C.  (1998).
   Algorithms for model-based Gaussian hierarchical clustering.
   \emph{SIAM Journal on Scientific Computing} 20:270-281.

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
   \link[mclust]{hc}
 }


 \examples{
 ## Load  Swiss banknotes data
 data(banknote)
 x <- banknote[,-1]

 ## Initial clusters with default arguments
 init <- InitClust(data = x, G = 2)
 print(init)

 ## Perform otrimle
 a <- otrimle(data = x, G = 2, initial = init,
              logicd = c(-Inf, -50, -10), ncores = 1)
 plot(a, what="clustering", data=x)
 }

